# === Multi-state intervals builder (onset, remission, relapse) ===
# Output: final_20_long_intervals_merged.parquet
# Schema preserved: person_id, start, stop, event, transition, age, sex

# ---------------- CONFIG ----------------
WORKSPACE_BUCKET = "xx"
OUT_DIR          = f"{WORKSPACE_BUCKET}/factorsobservation/survival_intervals_build"
OUT_PATH         = f"{OUT_DIR}/final_20_long_intervals_merged.parquet"

DEMO_URI   = f"{WORKSPACE_BUCKET}/factorsobservation/inputs/person.parquet"
RX_URI     = f"{WORKSPACE_BUCKET}/factorsobservation/inputs/drug.parquet"
CONDS_URI  = f"{WORKSPACE_BUCKET}/factorsobservation/inputs/condition_with_externalizing.parquet"
FEATURES_URI = f"{WORKSPACE_BUCKET}/factorsobservation/final_30_features_merged.parquet"

ORIGIN_MODE  = "features_min_start"   # "features_min_start" | "first_claim" | "any_rx"
TV_FREQ      = "14D"
MAX_INTERVALS = 120
CENSOR_DATE  = "2025-10-20"           # UTC
RELAPSE_GAP_DAYS        = 14
REMISSION_LOOKBACK_DAYS = 0

# ---------------- SPEED KNOBS ----------------
# Use all CPU threads and keep Arrow streaming/vectorized
import os, time, warnings, sys, subprocess
warnings.filterwarnings("ignore", category=FutureWarning)

os.environ.setdefault("POLARS_MAX_THREADS", "32")
os.environ.setdefault("ARROW_NUM_THREADS",  "32")
os.environ.setdefault("AWS_EC2_METADATA_DISABLED", "true")  # avoids fsspec stalls in some envs

# Lightweight installs (safe to re-run)
def _pip(p):
    try: __import__(p.split("==")[0].split(">=")[0].replace("-","_"))
    except Exception: subprocess.check_call([sys.executable, "-m", "pip", "install", p, "-q"])
for p in ("gcsfs>=2024.6.0", "pyarrow>=14.0.0", "polars>=1.6.0", "pandas>=2.2.2"):
    _pip(p)

import pandas as pd, numpy as np, polars as pl
storage_opts = {"token": "cloud"}

# ---------------- Utilities ----------------
def ensure_utc(ts):
    return pd.to_datetime(ts, errors="coerce", utc=True)

def _icd10_core(x: str) -> str:
    if not isinstance(x, str): return ""
    x = x.strip().upper()
    out = []
    for ch in x:
        if (ch == 'F' and not out) or ch.isdigit() or ch == '.':
            out.append(ch)
        else:
            break
    return "".join(out)

# OUD/Remission code sets 
OUD_SNOMED_CODES = {75544000, 5602001, 724653003, 191819002, 191820008}
OUD_SNOMED_FALLBACK = {"31948003","431855005","304521000000109","364146000","307824008","28947003","307825009","324121000000109"}
REMISSION_CODES_ICD10 = {"F11.11","F11.21","F11.91"}

# ---------------- Fast readers (project only needed columns) ----------------
def read_demo(uri=DEMO_URI):
    # person_id, date_of_birth, sex_at_birth_concept_id (optional)
    lf = pl.scan_parquet(uri).select([
        pl.col("person_id"),
        pl.col("date_of_birth"),
        pl.col("sex_at_birth_concept_id").alias("sex_at_birth_concept_id")
    ])
    df = lf.collect(engine="streaming").to_pandas()
    df["date_of_birth"] = ensure_utc(df["date_of_birth"])
    # sex: -1 unknown, 0 male(8507), 1 female(8532)
    if "sex_at_birth_concept_id" in df.columns:
        df["sex"] = df["sex_at_birth_concept_id"].map({8507:0, 8532:1}).fillna(-1).astype("int8")
    else:
        df["sex"] = -1
    return df[["person_id","date_of_birth","sex"]]

def read_rx(uri=RX_URI):
    # person_id, drug_exposure_start_datetime
    lf = pl.scan_parquet(uri).select(["person_id","drug_exposure_start_datetime"])
    df = lf.collect(engine="streaming").to_pandas()
    df["drug_exposure_start_datetime"] = ensure_utc(df["drug_exposure_start_datetime"])
    return df.dropna(subset=["drug_exposure_start_datetime"])

def read_conds(uri=CONDS_URI):
    # Only columns we touch
    lf = pl.scan_parquet(uri).select([
        "person_id",
        "condition_start_datetime",
        pl.col("standard_concept_name").cast(pl.Utf8, strict=False),
        pl.col("standard_concept_code").cast(pl.Utf8, strict=False),
        pl.col("condition_source_value").cast(pl.Utf8, strict=False),
        pl.col("source_concept_code").cast(pl.Utf8, strict=False),
    ])
    df = lf.collect(engine="streaming").to_pandas()
    df["condition_start_datetime"] = ensure_utc(df["condition_start_datetime"])
    return df.dropna(subset=["condition_start_datetime"])

# ---------------- OUD classification ----------------
def classify_oud_rows(conds: pd.DataFrame) -> pd.DataFrame:
    c = conds.copy()
    for col in ("standard_concept_name","standard_concept_code","condition_source_value","source_concept_code"):
        if col not in c.columns: c[col] = ""
    c["standard_concept_name"] = c["standard_concept_name"].fillna("").astype(str).str.lower()
    c["standard_concept_code"] = c["standard_concept_code"].fillna("").astype(str)
    c["source_concept_code"]   = c["source_concept_code"].fillna("").astype(str)
    c["icd10_core"]            = c["condition_source_value"].fillna("").astype(str).map(_icd10_core)

    is_icd_rem  = c["icd10_core"].isin(REMISSION_CODES_ICD10)
    is_name_rem = c["standard_concept_name"].str.contains("in remission", na=False)
    is_sno_rem  = c["standard_concept_code"].isin({str(x) for x in OUD_SNOMED_CODES})
    c["is_remission"] = (is_icd_rem | is_name_rem | is_sno_rem)

    oud_name = c["standard_concept_name"].str.contains(
        r"(?:\bopioid\b.*(?:use|depend|abuse|intoxication|withdrawal|overdose))", regex=True, na=False
    )
    oud_icd  = c["icd10_core"].str.startswith("F11") & ~c["is_remission"]
    oud_sno  = c["standard_concept_code"].isin({str(x) for x in OUD_SNOMED_CODES}) \
               | c["source_concept_code"].isin(OUD_SNOMED_FALLBACK)
    c["is_oud_positive"] = (oud_name | oud_icd | oud_sno) & ~c["is_remission"]

    out = c[["person_id","condition_start_datetime","is_oud_positive","is_remission"]].rename(
        columns={"condition_start_datetime":"dt"}
    )
    return out

def compute_oud_timeline(conds: pd.DataFrame):
    lab = classify_oud_rows(conds).sort_values(["person_id","dt"])
    first_oud = (lab.loc[lab.is_oud_positive, ["person_id","dt"]]
                   .groupby("person_id", as_index=False).min()
                   .rename(columns={"dt":"first_oud_date"}))
    # remission only after first OUD
    rem = lab.loc[lab.is_remission].merge(first_oud, on="person_id", how="left")
    rem = rem.loc[rem["dt"] > rem["first_oud_date"]]
    if REMISSION_LOOKBACK_DAYS > 0:
        win = pd.Timedelta(days=REMISSION_LOOKBACK_DAYS)
        pos = lab.loc[lab.is_oud_positive, ["person_id","dt"]]
        rem = rem.merge(pos, on="person_id", how="left", suffixes=("","_prev"))
        rem = rem[(rem["dt_prev"].isna()) | (rem["dt_prev"] <= rem["dt"] - win)]
        rem = rem.drop(columns=["dt_prev"])
    first_rem_after_oud = (rem[["person_id","dt"]]
                             .groupby("person_id", as_index=False).min()
                             .rename(columns={"dt":"first_remission_date"}))
    # relapse after remission + gap
    gap = pd.Timedelta(days=RELAPSE_GAP_DAYS)
    pos = lab.loc[lab.is_oud_positive, ["person_id","dt"]]
    rr  = pos.merge(first_rem_after_oud, on="person_id", how="inner")
    rr  = rr.loc[rr["dt"] >= rr["first_remission_date"] + gap]
    first_relapse_after_rem = (rr[["person_id","dt"]]
                                 .groupby("person_id", as_index=False).min()
                                 .rename(columns={"dt":"first_relapse_date"}))
    return first_oud, first_rem_after_oud, first_relapse_after_rem

# ---------------- Baseline/origin ----------------
def baseline_from_features(features_uri: str):
    lf = (pl.scan_parquet(features_uri)
            .select(["person_id","start"])
            .with_columns(pl.col("start").cast(pl.Datetime(time_unit="us", time_zone="UTC")))
            .group_by("person_id")
            .agg(pl.col("start").min().alias("origin")))
    return lf.collect(engine="streaming").to_pandas()


def compute_baseline(demo: pd.DataFrame, rx: pd.DataFrame, conds: pd.DataFrame,
                     origin_mode: str, features_uri: str|None) -> pd.DataFrame:
    frames = []
    mode = origin_mode.lower().strip()

    if mode in ("any_rx","first_claim"):
        if not rx.empty:
            first_rx = (rx.groupby("person_id", as_index=False)["drug_exposure_start_datetime"]
                          .min().rename(columns={"drug_exposure_start_datetime":"candidate"}))
            frames.append(first_rx)

    if mode == "features_min_start":
        if not features_uri:
            raise RuntimeError("ORIGIN_MODE=features_min_start requires FEATURES_URI")
        feat_min = baseline_from_features(features_uri).rename(columns={"origin":"candidate"})
        frames.append(feat_min)

    if mode in ("first_claim",):
        if not conds.empty:
            first_cond = (conds.groupby("person_id", as_index=False)["condition_start_datetime"]
                            .min().rename(columns={"condition_start_datetime":"candidate"}))
            frames.append(first_cond)

    if not frames:
        raise RuntimeError(f"compute_baseline: no candidates for ORIGIN_MODE={mode}")

    base = pd.concat(frames, ignore_index=True)
    base = (base.groupby("person_id", as_index=False)["candidate"].min()
                 .rename(columns={"candidate":"origin"}))
    # keep persons present in demo
    if demo is not None and not demo.empty:
        base = base.merge(demo[["person_id"]].drop_duplicates(), on="person_id", how="inner")
    return base

# ---------------- Multi-state base ----------------
def compute_multi_state_base(demo: pd.DataFrame, rx: pd.DataFrame, conds: pd.DataFrame) -> pd.DataFrame:
    base0 = compute_baseline(demo, rx, conds, ORIGIN_MODE, FEATURES_URI)  # person_id, origin

    first_oud, first_rem, first_rel = compute_oud_timeline(conds)

    # 0→1 onset
    b01 = (base0
           .merge(first_oud, on="person_id", how="left")
           .merge(demo, on="person_id", how="left"))
    b01["event_date"] = b01["first_oud_date"]
    b01["transition"] = "oud_onset"

    # 1→2 remission
    b12 = (first_oud
           .merge(first_rem, on="person_id", how="left")
           .merge(demo, on="person_id", how="left"))
    b12 = b12.dropna(subset=["first_oud_date"]).rename(columns={"first_oud_date":"origin"})
    b12["event_date"] = b12["first_remission_date"]
    b12["transition"] = "to_remission"

    # 2→3 relapse
    b23 = (first_rem
           .merge(first_rel, on="person_id", how="left")
           .merge(demo, on="person_id", how="left"))
    b23 = b23.dropna(subset=["first_remission_date"]).rename(columns={"first_remission_date":"origin"})
    b23["event_date"] = b23["first_relapse_date"]
    b23["transition"] = "to_relapse"

    base = pd.concat([b01, b12, b23], ignore_index=True, sort=False)

    base["event"]    = base["event_date"].notna().astype("int8")
    base["end_date"] = ensure_utc(base["event_date"].fillna(CENSOR_DATE))
    base["origin"]   = ensure_utc(base["origin"])
    base = base.loc[base["end_date"] >= base["origin"]].copy()

    base["age0"] = ((base["origin"] - base["date_of_birth"]).dt.days / 365.25).astype("float32")

    return base[["person_id","origin","end_date","event","age0","sex","transition"]]

# ---------------- Interval builder — Polars (vectorized, no Python loops) ----------------
def make_intervals_fast(base_df: pd.DataFrame, tv_freq: str, max_intervals: int) -> pd.DataFrame:
    b = pl.from_pandas(base_df)
    b = b.with_columns([
        pl.col("origin").cast(pl.Datetime(time_unit="us", time_zone="UTC")),
        pl.col("end_date").cast(pl.Datetime(time_unit="us", time_zone="UTC")),
        pl.col("event").cast(pl.Int8),
        pl.col("sex").cast(pl.Int8),
        pl.col("age0").cast(pl.Float32),
        pl.col("transition").cast(pl.Utf8),
    ])

    # Build cut grid per row (includes end bound)
    cuts = (
        b.with_columns(
            pl.datetime_ranges(
                start=pl.col("origin"),
                end=pl.col("end_date"),
                interval=tv_freq,
                closed="both"
            ).alias("cut")
        )
        .explode("cut")
        .with_columns([
            # previous cut within (person_id, transition); first one falls back to origin
            pl.coalesce([
                pl.col("cut").shift(1).over(["person_id", "transition"]),
                pl.col("origin")
            ]).alias("start"),
            pl.col("cut").alias("stop"),
        ])
        .drop("cut")
        .filter(pl.col("stop") > pl.col("start"))
    )

    # Cap number of intervals per (person_id, transition)
    # Robust per-(person_id, transition) running index 0,1,2,...
    if hasattr(pl, "row_number"):
        idx_expr = pl.row_number().over(["person_id","transition"])
    else:
    # Older Polars: cum_count needs a column name; any non-null column works
        idx_expr = pl.cum_count("start").over(["person_id","transition"])

    cuts = (
        cuts
        .with_columns(idx_expr.alias("__idx"))
        .filter(pl.col("__idx") < max_intervals)
        .drop("__idx")
    )


    # Put event only on the last interval per (person_id, transition)
    cuts = cuts.with_columns([
        (pl.col("stop") == pl.col("stop").max().over(["person_id","transition"]))
        .cast(pl.Int8).alias("__is_last")
    ])
    cuts = cuts.with_columns([
        pl.when(pl.col("__is_last") == 1).then(pl.col("event")).otherwise(0)
        .cast(pl.Int8).alias("event")
    ]).drop("__is_last")

    # Age at interval start: age0 + elapsed years since origin (timestamps are µs)
    cuts = cuts.with_columns([
        (
            pl.col("age0") +
            (pl.col("start").cast(pl.Int64) - pl.col("origin").cast(pl.Int64)) /
            pl.lit(365.25*24*3600*1_000_000)
        ).cast(pl.Float32).alias("age")
    ])

    out = cuts.select([
        "person_id",
        pl.col("start").alias("start"),
        pl.col("stop").alias("stop"),
        pl.col("event").alias("event"),
        pl.col("transition").alias("transition"),
        pl.col("age").alias("age"),
        pl.col("sex").alias("sex"),
    ]).sort(["person_id","transition","start"])

    return out.to_pandas()



# ---------------- RUN ----------------
t0 = time.time()
print("[1/5] Reading inputs (projection pushdown)…")
demo  = read_demo(DEMO_URI)
rx    = read_rx(RX_URI)
conds = read_conds(CONDS_URI)

print(f"[2/5] Building base (ORIGIN_MODE={ORIGIN_MODE})…")
base = compute_multi_state_base(demo, rx, conds)
print(f"    base rows: {len(base):,}")

print("[3/5] Making intervals (vectorized Polars)…")
long = make_intervals_fast(base, TV_FREQ, MAX_INTERVALS)
print(f"    intervals rows: {len(long):,}")

print("[4/5] Writing to:", OUT_PATH)
# Write via pandas (fsspec to GCS). Keeps same file format as before.
long.to_parquet(OUT_PATH, index=False, storage_options=storage_opts)

print("[5/5] Sanity checks (Polars scan)…")
d = pl.scan_parquet(OUT_PATH)
print(
    d.select([
        pl.len().alias("rows"),
        pl.col("person_id").n_unique().alias("unique_ids"),
        pl.col("event").sum().alias("events_total"),
        pl.col("start").min().alias("min_start"),
        pl.col("stop").max().alias("max_stop"),
    ]).collect(streaming=True)
)
print(
    d.group_by("transition").agg([
        pl.len().alias("rows"),
        pl.col("person_id").n_unique().alias("unique_ids"),
        pl.col("event").sum().alias("events_total"),
    ]).sort("events_total", descending=True).collect(streaming=True)
)

print(f"Done in {time.time()-t0:,.1f}s.")
