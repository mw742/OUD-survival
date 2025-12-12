#!/usr/bin/env python3
"""
OUD → Gene → Drug ranking with:
  • pathway-bridging,
  • alias-aware labeling (parent/child salt & brand variants),
  • per-drug "why" provenance,
  • deduplicated labels.
Inputs
------
nodes.csv:
  id,node_type,name,[approved],[bbb],[research],...
edges.csv:
  src,dst,edge_type,src_type,dst_type,[pchembl],[moa],[phase],
  [paper1_gws],[paper2_gws],...
  [replication_count],[posterior],[p_rank],...

Key behaviors
-------------
• ASSOCIATED_WITH Gene→OUD becomes HAS_GENE OUD→Gene with evidence weight.
• TARGETS (any orientation) normalized to TARGETED_BY Gene→Drug with potency/flags weight.
• IN_PATHWAY Gene→Pathway used to add OUD→Gene HAS_GENE_BRIDGED edges (γ penalty).
• Pathway provenance recorded: which pathways bridged a given gene.
• Alias-aware labels:
    alias-mode borrow (default): if a child drug has no targets, borrow labels from its parent.
    alias-mode fanout          : copy parent’s Gene→Drug edges to all children present in nodes.
    alias-mode none            : no alias help.
• Output files:
    --out <TSV>              : main rankings with `gene_labels` and `why` columns.
    --why-out <TSV> (opt)    : long-form provenance per (drug, gene).
"""

import argparse, csv, os, sys, math, logging
from collections import defaultdict, Counter

# ----- graph-tool -----
try:
    from graph_tool.all import Graph, GraphView, pagerank
except Exception as e:
    print("[FATAL] Need graph-tool in your environment.\n"
          f"Import error: {e}", file=sys.stderr)
    sys.exit(1)

# Utilities

def fail(msg: str) -> None:
    logging.error(msg); sys.exit(2)

def str_to_bool(x) -> bool:
    if x is None: return False
    s = str(x).strip().lower()
    return s in {"1","true","t","yes","y"}

def coerce_float(x, default=None):
    try:
        if x is None or str(x).strip()=="":
            return default
        return float(x)
    except Exception:
        return default

def coerce_int(x, default=None):
    try:
        if x is None or str(x).strip()=="":
            return default
        return int(x)
    except Exception:
        return default

def read_csv_rows(path: str) -> list[dict]:
    if not os.path.exists(path):
        fail(f"CSV not found: {path}  (cwd={os.getcwd()})")
    with open(path, newline="", encoding="utf-8-sig") as f:
        rdr = csv.DictReader(f)
        return [{(k or "").strip(): (v.strip() if isinstance(v,str) else v) for k,v in row.items()} for row in rdr]

def infer_type_from_id(nid: str) -> str:
    if nid.upper() == "OUD": return "Disease"
    u = nid.upper()
    if u.startswith("CHEMBL") or u in {"NS309","1-EBIO"}: 
        return "Drug"
    if u.startswith("R-HSA-") or u.startswith("REACTOME:"):
        return "Pathway"
    return "Gene"

# Meta (alias) helpers

def load_meta_alias(meta_csv: str):
    """
    Returns:
      child2parent: dict child_id -> parent_id
      parent_name:  dict parent_id -> parent_name (fallback to id if missing)
      cluster:      dict parent_id -> set of all ids in that parent cluster (incl parent)
    Accepts columns: (child_id,parent_id,parent_name?) or (drug_id,parent_id,parent_name?)
    """
    if not meta_csv:
        return {}, {}, {}
    rows = read_csv_rows(meta_csv)
    child2parent = {}
    parent_name = {}
    cluster = defaultdict(set)
    for r in rows:
        child = r.get("child_id") or r.get("drug_id")
        parent = r.get("parent_id") or child
        if not child: 
            continue
        child = child.strip()
        parent = (parent or child).strip()
        child2parent[child] = parent
        if r.get("parent_name"): parent_name[parent] = r["parent_name"]
        cluster[parent].add(child)
        cluster[parent].add(parent)
    return child2parent, parent_name, cluster

def parent_of(drug_id: str, child2parent: dict) -> str:
    return child2parent.get(drug_id, drug_id)

# Edge weighting

def oud_assoc_weight(row: dict) -> float:
    repl = coerce_int(row.get("replication_count"))
    if repl is not None:
        return max(1e-9, min(1.5, 0.5 + 0.3 * repl))

    if str_to_bool(row.get("paper1_mvp_oud_gws")): return 1.0
    if str_to_bool(row.get("paper2_oud_only_gws")): return 0.8
    if str_to_bool(row.get("paper2_mtag_gws")):     return 0.5

    post = coerce_float(row.get("posterior"))
    if post is not None:
        return max(1e-9, min(1.2, 0.2 + post))

    p_rank = coerce_float(row.get("p_rank"))
    if p_rank is not None:
        return max(1e-9, 1.0 / math.log(3.0 + p_rank))

    return 1e-9

def target_edge_weight(drug_v, v_approved, v_bbb, v_research, row: dict) -> float:
    w = 1.0
    pchembl = coerce_float(row.get("pchembl"))
    if pchembl is not None:
        w += max(0.0, min(1.0, (pchembl - 5.0) / 4.0))
    moa = (row.get("moa") or "").lower()
    if any(k in moa for k in ("antagonist","inverse","inhibitor")):
        w += 0.2
    if any(k in moa for k in ("agonist","partial agonist","modulator")):
        w += 0.1
    phase = row.get("phase")
    if str(phase).isdigit() and int(phase) >= 3:
        w += 0.3
    if v_approved[drug_v]: w += 0.7
    if v_bbb[drug_v]:      w += 0.2
    if v_research[drug_v]: w -= 0.5
    return max(w, 1e-9)

# Graph building

def build_graph(nodes_csv: str, edges_csv: str, gamma: float, max_bridged: int, min_seeds_per_pathway: int,
                alias_mode: str, meta_csv: str):
    # alias meta
    child2parent, parent_name, cluster = load_meta_alias(meta_csv)

    nodes = read_csv_rows(nodes_csv)
    edges = read_csv_rows(edges_csv)

    g = Graph(directed=True)
    # vertex props
    v_id   = g.new_vertex_property("string")
    v_type = g.new_vertex_property("string")
    v_name = g.new_vertex_property("string")
    v_approved = g.new_vertex_property("bool")
    v_bbb      = g.new_vertex_property("bool")
    v_research = g.new_vertex_property("bool")
    # attach early
    g.vp["id"], g.vp["node_type"], g.vp["name"] = v_id, v_type, v_name
    g.vp["approved"], g.vp["bbb"], g.vp["research"] = v_approved, v_bbb, v_research

    idx: dict[str, "Vertex"] = {}

    def touch_v(nid: str, ntype: str="", name: str=""):
        v = idx.get(nid)
        if v is None:
            v = g.add_vertex()
            idx[nid] = v
            v_id[v] = nid
            v_type[v] = ntype or ""
            v_name[v] = name or nid
            # defaults
            name_l = v_name[v].lower()
            is_mat = name_l in {"naltrexone","methadone","buprenorphine"}
            v_approved[v] = is_mat
            v_bbb[v]      = False
            v_research[v] = ("research" in name_l)
        else:
            if ntype and not v_type[v]:
                v_type[v] = ntype
            if name and (not v_name[v] or v_name[v] == v_id[v]):
                v_name[v] = name
        return v

    # create nodes
    for n in nodes:
        nid  = n.get("id") or n.get("node_id") or ""
        if not nid: 
            logging.warning("Skipping node with empty id"); 
            continue
        ntype = n.get("node_type") or ""
        name  = n.get("name") or parent_name.get(nid) or nid
        v = touch_v(nid, ntype, name)
        if "approved" in n: v_approved[v] = str_to_bool(n["approved"])
        if "bbb" in n:      v_bbb[v]      = str_to_bool(n["bbb"])
        if "research" in n: v_research[v] = str_to_bool(n["research"])

    # edge props
    e_type = g.new_edge_property("string")
    e_w    = g.new_edge_property("double", val=1.0)
    g.ep["edge_type"], g.ep["weight"] = e_type, e_w

    # pathway membership & provenance maps
    path_members: dict[str, set[str]] = defaultdict(set)   # pathway_id -> {gene_id}
    path_names:   dict[str, str] = {}                      # pathway_id -> name
    bridged_by_paths: dict[int, set[str]] = defaultdict(set)  # gene_idx -> set(pathway_id)

    # ingest edges
    for e in edges:
        s_id, t_id = (e.get("src","") or "").strip(), (e.get("dst","") or "").strip()
        if not s_id or not t_id: 
            continue
        et = (e.get("edge_type") or "").strip().upper()
        s_type = (e.get("src_type") or "").strip() or infer_type_from_id(s_id)
        t_type = (e.get("dst_type") or "").strip() or infer_type_from_id(t_id)
        if t_id.upper() == "OUD": t_type = "Disease"
        if s_id.upper() == "OUD": s_type = "Disease"

        vs = touch_v(s_id, s_type, "")
        vt = touch_v(t_id, t_type, "")

        # raw edge (for provenance)
        e1 = g.add_edge(vs, vt)
        e_type[e1] = et
        e_w[e1]    = 1.0

        # pathway membership (assume Gene -> Pathway)
        if et == "IN_PATHWAY":
            path_members[t_id].add(s_id)
            if t_id not in path_names:
                path_names[t_id] = v_name[vt] or t_id
            continue

        # Gene→OUD association
        if et == "ASSOCIATED_WITH" and (t_id.upper()=="OUD" or v_type[vt]=="Disease"):
            w = oud_assoc_weight(e)
            e2 = g.add_edge(vt, vs)  # OUD -> Gene
            e_type[e2] = "HAS_GENE"
            e_w[e2]    = max(w, 1e-9)
            continue

        # Normalize TARGETS
        if et in {"TARGETS","TARGETED_BY","GENE_TARGET"}:
            # orientation
            if v_type[vs]=="Drug" and v_type[vt]=="Gene":
                gene_v, drug_id = vt, v_id[vs]
            elif v_type[vs]=="Gene" and v_type[vt]=="Drug":
                gene_v, drug_id = vs, v_id[vt]
            else:
                # last try infer
                if infer_type_from_id(s_id)=="Drug" and infer_type_from_id(t_id)=="Gene":
                    gene_v, drug_id = vt, s_id
                elif infer_type_from_id(s_id)=="Gene" and infer_type_from_id(t_id)=="Drug":
                    gene_v, drug_id = vs, t_id
                else:
                    continue

            # which drug(s) to attach: alias fanout or single
            attach_ids = [drug_id]
            if alias_mode == "fanout":
                root = parent_of(drug_id, child2parent)
                if root in cluster:
                    # only attach to drugs that exist in the graph (present in nodes/edges)
                    attach_ids = [d for d in cluster[root] if d in idx]
                else:
                    attach_ids = [drug_id]

            for did in attach_ids:
                dv = touch_v(did, "Drug", parent_name.get(did, ""))
                w = target_edge_weight(dv, v_approved, v_bbb, v_research, e)
                e2 = g.add_edge(gene_v, dv)
                e_type[e2] = "TARGETED_BY"
                e_w[e2]    = w
            continue

    # find OUD vertex
    OUD_V = None
    for v in g.vertices():
        if v_id[v].upper() == "OUD":
            OUD_V = v; break
    if OUD_V is None:
        fail("No node with id='OUD' found.")

    # collect direct OUD gene weights
    direct_oud_gene_w: dict[int, float] = {}
    for e in OUD_V.out_edges():
        if e_type[e] == "HAS_GENE":
            direct_oud_gene_w[int(e.target())] = e_w[e]

    # bridge via pathways
    bridged_count = 0
    for pw_id, members in path_members.items():
        direct_members = []
        for gsym in members:
            gv = idx.get(gsym)
            if gv is not None and int(gv) in direct_oud_gene_w:
                direct_members.append(gv)
        if len(direct_members) < min_seeds_per_pathway:
            continue
        agg = sum(direct_oud_gene_w[int(v)] for v in direct_members)
        base_bonus = max(1e-9, gamma * (agg / max(1, len(members))))
        for gsym in members:
            if bridged_count >= max_bridged: break
            gv = idx.get(gsym)
            if gv is None or gv == OUD_V: continue
            if int(gv) in direct_oud_gene_w: continue  # already direct
            e2 = g.add_edge(OUD_V, gv)
            e_type[e2] = "HAS_GENE_BRIDGED"
            e_w[e2]    = base_bonus
            bridged_count += 1
            # record pathway provenance
            bridged_by_paths[int(gv)].add(pw_id)

    logging.info(f"Pathway-bridged genes added: {bridged_count}")

    # keep handy maps on g (for later access through GraphView)
    g.graph_properties["_bridged_by_paths"] = g.new_graph_property("object")
    g.graph_properties["_bridged_by_paths"] = bridged_by_paths
    g.graph_properties["_path_names"] = g.new_graph_property("object")
    g.graph_properties["_path_names"] = path_names
    g.graph_properties["_child2parent"] = g.new_graph_property("object")
    g.graph_properties["_child2parent"] = child2parent
    g.graph_properties["_cluster"] = g.new_graph_property("object")
    g.graph_properties["_cluster"] = cluster

    return g

# Subgraph + normalization

def project_dgd(g: Graph) -> GraphView:
    keep_v = g.new_vertex_property("bool")
    for v in g.vertices():
        keep_v[v] = g.vp["node_type"][v] in {"Disease","Gene","Drug"}
    keep_e = g.new_edge_property("bool")
    for e in g.edges():
        keep_e[e] = g.ep["edge_type"][e] in {"HAS_GENE","HAS_GENE_BRIDGED","TARGETED_BY"}
    return GraphView(g, vfilt=keep_v, efilt=keep_e, directed=True)

def normalize_outgoing_weights(g: GraphView) -> None:
    w = g.ep["weight"]
    out_sum = g.new_vertex_property("double", val=0.0)
    for e in g.edges(): out_sum[e.source()] += w[e]
    for e in g.edges():
        denom = out_sum[e.source()]
        if denom > 0: w[e] = w[e] / denom

# Personalization

def build_personalization(sg: GraphView, source_id="OUD", approved_prior=0.5):
    pers = sg.new_vertex_property("double", val=0.0)
    oud = None
    for v in sg.vertices():
        if sg.vp["id"][v] == source_id:
            oud = v; break
    if oud is None: fail(f"Source node {source_id} not found.")
    pers[oud] = 1.0
    for e in oud.out_edges():
        if sg.ep["edge_type"][e] in {"HAS_GENE","HAS_GENE_BRIDGED"}:
            pers[e.target()] += sg.ep["weight"][e]
    for v in sg.vertices():
        if sg.vp["node_type"][v] == "Drug" and sg.vp["approved"][v]:
            pers[v] += approved_prior
    total = sum(pers[v] for v in sg.vertices())
    if total <= 0: fail("Personalization vector has zero mass.")
    for v in sg.vertices(): pers[v] = pers[v] / total
    return pers

# Ranking + explanations

def _oud_gene_maps(sg: GraphView):
    """Return maps: status{gene_idx:'direct'|'bridged'}, oud_w{gene_idx:weight}."""
    status = {}
    oud_w  = {}
    oud = next(v for v in sg.vertices() if sg.vp["id"][v] == "OUD")
    for e in oud.out_edges():
        et = sg.ep["edge_type"][e]
        if et in {"HAS_GENE","HAS_GENE_BRIDGED"}:
            gidx = int(e.target())
            status[gidx] = "direct" if et=="HAS_GENE" else "bridged"
            oud_w[gidx]  = float(sg.ep["weight"][e])
    return status, oud_w

def _best_status_merge(cur: str|None, new: str) -> str:
    order = {"direct":3, "bridged":2, "nonseed":1, None:0}
    return new if order[new] > order.get(cur,0) else cur

def rank_with_why(sg: GraphView, pr, alias_mode: str, borrow_from_parent: bool, why_top: int, why_out_path: str|None):
    """
    Returns rows:
      (drug_name, drug_id, score, oud_genes[], approved, bbb, gene_labels[], why_str)
    Also writes a long-form provenance TSV if why_out_path is provided.
    """
    bridged_by_paths = sg.base.graph_properties["_bridged_by_paths"]
    path_names       = sg.base.graph_properties["_path_names"]
    child2parent     = sg.base.graph_properties["_child2parent"]
    cluster          = sg.base.graph_properties["_cluster"]

    status_map, oud_w_map = _oud_gene_maps(sg)

    # Per-drug provenance collection (for why_out)
    long_rows = []

    rows = []
    w_prop = sg.ep["weight"]

    # Build fast map: for each drug, collect inbound Gene edges
    for v in sg.vertices():
        if sg.vp["node_type"][v] != "Drug":
            continue
        did = sg.vp["id"][v]
        dname = sg.vp["name"][v]
        approved = bool(sg.vp["approved"][v]); bbb = bool(sg.vp["bbb"][v])

        # Gather Gene→Drug edges
        in_edges = []
        for e in v.in_edges():
            if sg.ep["edge_type"][e] == "TARGETED_BY":
                gidx = int(e.source())
                gsym = sg.vp["id"][e.source()]
                tw   = float(w_prop[e])
                in_edges.append((gidx, gsym, tw))

        # If none, try "borrow" from parent for labeling/explanations
        borrowed = False
        borrowed_edges = []
        if not in_edges and borrow_from_parent:
            root = child2parent.get(did, did)
            # pick any cluster member present that *has* targets (prefer the root)
            candidates = []
            if root in cluster:
                for mem in cluster[root]:
                    # find that vertex if present
                    for vv in sg.vertices():
                        if sg.vp["id"][vv] == mem and sg.vp["node_type"][vv]=="Drug":
                            # collect its inbound gene edges
                            tmp = []
                            for e in vv.in_edges():
                                if sg.ep["edge_type"][e]=="TARGETED_BY":
                                    gidx = int(e.source())
                                    gsym = sg.vp["id"][e.source()]
                                    tw   = float(w_prop[e])
                                    tmp.append((gidx, gsym, tw))
                            if tmp:
                                candidates.append((vv, tmp))
            if candidates:
                # choose the candidate with most inbound edges
                vv, tmp = max(candidates, key=lambda x: len(x[1]))
                in_edges = tmp
                borrowed = True

        # Collapse duplicates by gene: merge multiple sources
        by_gene = {}
        for gidx, gsym, tw in in_edges:
            if gidx not in by_gene:
                by_gene[gidx] = {"sym": gsym, "tw": 0.0}
            by_gene[gidx]["tw"] += tw  # sum parallel weights

        # Split into OUD genes vs nonseed
        oud_genes_syms = []
        labels_map = {}  # gsym -> best_status
        reasons = []     # tuples for sorting: (contrib, text)

        for gidx, rec in by_gene.items():
            gsym = rec["sym"]; tw = rec["tw"]
            st  = status_map.get(gidx)
            # status & OUD→Gene weight
            if st:
                oud_w = oud_w_map.get(gidx, 0.0)
                oud_genes_syms.append(gsym)
                labels_map[gsym] = _best_status_merge(labels_map.get(gsym), st)
                # provenance text
                if st == "direct":
                    why_txt = f"{gsym}(direct; w_oud={oud_w:.3f}, w_tgt={tw:.3f})"
                else:
                    # list a few pathways
                    pws = sorted(bridged_by_paths.get(gidx, []))
                    if pws:
                        pretty = [f"{pid}({path_names.get(pid, pid)})" for pid in pws[:3]]
                        extra  = " via " + "; ".join(pretty)
                    else:
                        extra = ""
                    why_txt = f"{gsym}(bridged{extra}; w_oud={oud_w:.3f}, w_tgt={tw:.3f})"
                contrib = (oud_w or 0.0) * tw
                reasons.append((contrib, why_txt))

                # append long-form row
                long_rows.append({
                    "drug_id": did, "drug_name": dname, "borrowed": int(borrowed),
                    "gene": gsym, "status": st,
                    "oud_w": f"{(oud_w or 0.0):.6f}",
                    "tgt_w": f"{tw:.6f}",
                    "contribution": f"{contrib:.6f}",
                    "pathways": ";".join(sorted(bridged_by_paths.get(gidx, [])))
                })
            else:
                # nonseed fallback
                labels_map[gsym] = _best_status_merge(labels_map.get(gsym), "nonseed")
                contrib = 0.0
                reasons.append((contrib, f"{gsym}(nonseed; w_tgt={tw:.3f})"))
                long_rows.append({
                    "drug_id": did, "drug_name": dname, "borrowed": int(borrowed),
                    "gene": gsym, "status": "nonseed",
                    "oud_w": f"{0.0:.6f}", "tgt_w": f"{tw:.6f}",
                    "contribution": f"{0.0:.6f}", "pathways": ""
                })

        # Build label list (deterministic order: direct>bridged>nonseed, then alpha)
        priority = {"direct":0, "bridged":1, "nonseed":2}
        label_items = sorted([(gsym, labels_map[gsym]) for gsym in labels_map.keys()],
                             key=lambda x: (priority[x[1]], x[0]))
        gene_labels = [f"{g}({s})" for g,s in label_items]

        # why string: top-N reasons by contribution, keep at least nonseed if nothing else
        reasons.sort(key=lambda x: x[0], reverse=True)
        why_list = [t for (_c, t) in reasons[:max(1, why_top)]]
        why_str = "; ".join(why_list) if why_list else ("approved_prior" if approved else "")

        rows.append((
            dname, did, float(pr[v]),
            sorted(set(oud_genes_syms)),
            approved, bbb,
            gene_labels, why_str
        ))

    rows.sort(key=lambda x: x[2], reverse=True)

    # write long-form provenance if requested
    if why_out_path:
        d = os.path.dirname(why_out_path)
        if d: os.makedirs(d, exist_ok=True)
        with open(why_out_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["drug_id","drug_name","borrowed","gene","status","oud_w","tgt_w","contribution","pathways"])
            for r in long_rows:
                w.writerow([r["drug_id"], r["drug_name"], r["borrowed"], r["gene"], r["status"],
                            r["oud_w"], r["tgt_w"], r["contribution"], r["pathways"]])

    return rows

def write_tsv(path: str, rows):
    if not path: return
    d = os.path.dirname(path)
    if d: os.makedirs(d, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "rank","drug_name","drug_id","score",
            "n_oud_genes","oud_genes","gene_labels","why","approved","bbb"
        ])
        for i,(name,did,score,oud_genes,approved,bbb,labels,why) in enumerate(rows, start=1):
            w.writerow([
                i, name, did, f"{score:.8f}",
                len(oud_genes),
                ",".join(oud_genes) if oud_genes else "",
                "; ".join(labels) if labels else "",
                why,
                int(approved), int(bbb)
            ])

# Main

def main():
    ap = argparse.ArgumentParser(description="Disease–Gene–Drug ranking with explanations and alias-aware labels")
    ap.add_argument("--nodes", required=True, help="nodes.csv")
    ap.add_argument("--edges", required=True, help="edges.csv")
    ap.add_argument("--meta", default="", help="Optional alias meta CSV (child_id,parent_id,...)")
    ap.add_argument("--alias-mode", default="borrow", choices=["borrow","fanout","none"],
                    help="How to use alias meta: borrow=use parent labels if no edges; fanout=copy edges to children; none=ignore")
    ap.add_argument("--alpha", type=float, default=0.85, help="PageRank damping")
    ap.add_argument("--gamma", type=float, default=0.20, help="Pathway-bridging penalty (0..1)")
    ap.add_argument("--max-bridged-genes", type=int, default=200, help="Cap for pathway-bridged OUD→Gene edges")
    ap.add_argument("--min-seeds-per-pathway", type=int, default=1, help="Minimum direct OUD genes required in a pathway to bridge others")
    ap.add_argument("--approved-prior", type=float, default=0.5, help="Personalization prior added per approved drug")
    ap.add_argument("--why-top", type=int, default=3, help="Top-N contributing genes to include in 'why' text")
    ap.add_argument("--why-out", default="", help="Optional TSV to dump long-form per-drug provenance")
    ap.add_argument("--topk", type=int, default=100, help="Top-K to print (console); TSV contains all")
    ap.add_argument("--out", default="", help="TSV path for rankings")
    ap.add_argument("--log", default="INFO", choices=["DEBUG","INFO","WARNING","ERROR"], help="Log level")
    args = ap.parse_args()

    logging.basicConfig(level=getattr(logging, args.log), format="[%(levelname)s] %(message)s")

    borrow_from_parent = (args.alias_mode == "borrow")
    g  = build_graph(args.nodes, args.edges, gamma=args.gamma,
                     max_bridged=args.max_bridged_genes,
                     min_seeds_per_pathway=args.min_seeds_per_pathway,
                     alias_mode=args.alias_mode, meta_csv=args.meta)

    sg = project_dgd(g)
    normalize_outgoing_weights(sg)

    pers = build_personalization(sg, source_id="OUD", approved_prior=args.approved_prior)
    pr   = pagerank(sg, damping=args.alpha, pers=pers, weight=sg.ep["weight"])

    rows = rank_with_why(sg, pr,
                         alias_mode=args.alias_mode,
                         borrow_from_parent=borrow_from_parent,
                         why_top=args.why_top,
                         why_out_path=(args.why_out or None))

    # console preview
    topk = rows[:max(1, args.topk)]
    print(f"All drugs ranked by weighted & normalized PPR (alpha={args.alpha}):")
    for name, did, score, oud_genes, approved, bbb, labels, why in topk:
        genes_str = ", ".join(oud_genes) if oud_genes else "-"
        flags = []
        if approved: flags.append("approved")
        if bbb:      flags.append("BBB")
        flag_str = f"  ({', '.join(flags)})" if flags else ""
        label_str = "; ".join(labels) if labels else "-"
        print(f"{name:<28} {did:<16} score={score:.6f}  oud_genes=[{genes_str}]  labels=[{label_str}]  why=[{why}]{flag_str}")

    write_tsv(args.out, rows)
    if args.out:
        print(f"\n[OK] Wrote: {args.out}")
    if args.why_out:
        print(f"[OK] Wrote provenance: {args.why_out}")

if __name__ == "__main__":
    main()

