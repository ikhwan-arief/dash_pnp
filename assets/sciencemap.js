/*! Copyright (c) 2026 Ikhwan Arief (ikhwan[at]eng.unand.ac.id). All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 * (https://polyformproject.org/licenses/noncommercial/1.0.0). Commercial use requires prior
 * written permission from the copyright holder. See LICENSE and NOTICE. */
// Citation maps in the manner of VOSviewer: association-strength normalisation, modularity clustering and the VOS layout.
// Everything is deterministic, so the same snapshot and settings always give the same map.
export const MAP_LIMITS = { nodes: 300, edges: 3000, commonReference: 500, cocitationPairs: 20000, minClusterSize: 3 };
export const MAP_KINDS = ['coupling', 'direct', 'cocitation'];
const pairKey = (a, b) => a < b ? `${a}|${b}` : `${b}|${a}`;
const add = (map, a, b, n = 1) => { if (a !== b) { const key = pairKey(a, b); map.set(key, (map.get(key) || 0) + n); } };

// Raw link counts between works. coupling: shared references; direct: one cites the other; cocitation: cited together.
// External co-citation pairs ([a, b, count], from all citing works) replace the in-snapshot count when available.
export function rawLinks(kind, works, ids, cocitationPairs = null) {
    const links = new Map();
    if (kind === 'coupling') {
        const citers = new Map();
        for (const w of works) if (ids.has(w.id))
            for (const ref of new Set(w.references || [])) { if (!citers.has(ref)) citers.set(ref, []); citers.get(ref).push(w.id); }
        // A reference shared by hundreds of works (a handbook, a standard) carries little structure and costs quadratic time.
        for (const list of citers.values()) if (list.length > 1 && list.length <= MAP_LIMITS.commonReference)
            for (let i = 0; i < list.length; i++) for (let j = i + 1; j < list.length; j++) add(links, list[i], list[j]);
    } else if (kind === 'direct') {
        for (const w of works) if (ids.has(w.id))
            for (const ref of new Set(w.references || [])) if (ids.has(ref)) add(links, w.id, ref);
    } else if (kind === 'cocitation') {
        if (cocitationPairs) { for (const [a, b, n] of cocitationPairs) if (ids.has(a) && ids.has(b)) add(links, a, b, n); }
        else for (const w of works) {
            const refs = [...new Set((w.references || []).filter(id => ids.has(id)))];
            for (let i = 0; i < refs.length; i++) for (let j = i + 1; j < refs.length; j++) add(links, refs[i], refs[j]);
        }
    } else throw new Error(`Unknown map type "${kind}".`);
    return links;
}

// Association strength s_ij = 2m·a_ij / (k_i·k_j), with k the total link strength and m the total of all links.
export function associationStrength(links) {
    const total = new Map(); let m = 0;
    for (const [key, n] of links) { const [a, b] = key.split('|'); total.set(a, (total.get(a) || 0) + n); total.set(b, (total.get(b) || 0) + n); m += n; }
    const out = new Map();
    for (const [key, n] of links) { const [a, b] = key.split('|'); out.set(key, 2 * m * n / (total.get(a) * total.get(b))); }
    return { strength: out, total };
}

// Louvain modularity clustering on a weighted graph given as n nodes and [i, j, w] edges; returns a cluster index per node.
export function louvain(n, edges, resolution = 1) {
    let community = Array.from({ length: n }, (_, i) => i), nodes = n, graph = edges.map(e => [...e]);
    const assignment = Array.from({ length: n }, (_, i) => i);
    for (let level = 0; level < 20; level++) {
        const adj = Array.from({ length: nodes }, () => new Map()), degree = new Float64Array(nodes); let m2 = 0;
        for (const [i, j, w] of graph) {
            adj[i].set(j, (adj[i].get(j) || 0) + w); if (i !== j) adj[j].set(i, (adj[j].get(i) || 0) + w);
            degree[i] += w; degree[j] += w; m2 += 2 * w;
        }
        if (!m2) break;
        community = Array.from({ length: nodes }, (_, i) => i);
        const tot = Float64Array.from(degree);
        let moved = true, any = false, rounds = 0;
        while (moved && rounds++ < 50) {
            moved = false;
            for (let i = 0; i < nodes; i++) {
                const own = community[i], weights = new Map();
                for (const [j, w] of adj[i]) if (j !== i) weights.set(community[j], (weights.get(community[j]) || 0) + w);
                tot[own] -= degree[i];
                let best = own, gain = (weights.get(own) || 0) - resolution * tot[own] * degree[i] / m2;
                for (const [c, w] of [...weights].sort((a, b) => a[0] - b[0])) {
                    const g = w - resolution * tot[c] * degree[i] / m2;
                    if (g > gain + 1e-12) { gain = g; best = c; }
                }
                tot[best] += degree[i];
                if (best !== own) { community[i] = best; moved = true; any = true; }
            }
        }
        const renumber = new Map(); community = community.map(c => { if (!renumber.has(c)) renumber.set(c, renumber.size); return renumber.get(c); });
        for (let v = 0; v < n; v++) assignment[v] = community[assignment[v]];
        if (!any || renumber.size === nodes) break;
        const merged = new Map();
        for (const [i, j, w] of graph) { const a = community[i], b = community[j], key = a <= b ? `${a}|${b}` : `${b}|${a}`; merged.set(key, (merged.get(key) || 0) + w); }
        graph = [...merged].map(([key, w]) => { const [a, b] = key.split('|').map(Number); return [a, b, w]; });
        nodes = renumber.size;
    }
    return assignment;
}

// Small deterministic generator for the starting layout.
function random(seed) { let s = seed >>> 0; return () => { s = (s + 0x6D2B79F5) >>> 0; let t = s; t = Math.imul(t ^ t >>> 15, t | 1); t ^= t + Math.imul(t ^ t >>> 7, t | 61); return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }

// VOS layout: minimise Σ s_ij·d_ij² − Σ d_ij over all pairs (attraction by similarity, uniform repulsion), then centre and scale to mean distance 1.
export function vosLayout(n, edges, { iterations = 600, seed = 1 } = {}) {
    const rnd = random(seed), x = new Float64Array(n), y = new Float64Array(n);
    for (let i = 0; i < n; i++) { x[i] = rnd() * 2 - 1; y[i] = rnd() * 2 - 1; }
    if (n < 2) return Array.from({ length: n }, () => [0, 0]);
    let step = 1;
    for (let it = 0; it < iterations; it++) {
        const gx = new Float64Array(n), gy = new Float64Array(n);
        for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
            const dx = x[i] - x[j], dy = y[i] - y[j], d = Math.hypot(dx, dy) || 1e-9;
            gx[i] -= dx / d; gy[i] -= dy / d; gx[j] += dx / d; gy[j] += dy / d;
        }
        for (const [i, j, s] of edges) { const dx = x[i] - x[j], dy = y[i] - y[j]; gx[i] += 2 * s * dx; gy[i] += 2 * s * dy; gx[j] -= 2 * s * dx; gy[j] -= 2 * s * dy; }
        let norm = 0; for (let i = 0; i < n; i++) norm = Math.max(norm, Math.hypot(gx[i], gy[i]));
        if (!norm) break;
        for (let i = 0; i < n; i++) { x[i] -= step * gx[i] / norm; y[i] -= step * gy[i] / norm; }
        step *= 0.985;
    }
    const cx = x.reduce((a, b) => a + b, 0) / n, cy = y.reduce((a, b) => a + b, 0) / n;
    let mean = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) mean += Math.hypot(x[i] - x[j], y[i] - y[j]);
    mean /= n * (n - 1) / 2 || 1;
    return Array.from({ length: n }, (_, i) => [(x[i] - cx) / mean, (y[i] - cy) / mean]);
}

// Full map: pick the most strongly linked works, normalise, cluster and lay out. label(work) and topicOf(work) come from the caller.
export function buildScienceMap(works, kind, { pool = works, minCitations = 0, maxNodes = MAP_LIMITS.nodes, cocitationPairs = null, label = w => w.id } = {}) {
    maxNodes = Math.max(10, Math.min(MAP_LIMITS.nodes, Number(maxNodes) || MAP_LIMITS.nodes));
    const candidates = pool.filter(w => (w.citations || 0) >= Number(minCitations || 0)), byId = new Map(candidates.map(w => [w.id, w]));
    const links = rawLinks(kind, works, new Set(byId.keys()), kind === 'cocitation' ? cocitationPairs : null);
    const tls = new Map();
    for (const [key, n] of links) for (const id of key.split('|')) tls.set(id, (tls.get(id) || 0) + n);
    const chosen = [...tls.keys()].sort((a, b) => tls.get(b) - tls.get(a) || (byId.get(b).citations || 0) - (byId.get(a).citations || 0) || (a < b ? -1 : 1)).slice(0, maxNodes);
    const keep = new Set(chosen), kept = new Map([...links].filter(([key]) => { const [a, b] = key.split('|'); return keep.has(a) && keep.has(b); }));
    // As VOSviewer's default, only the largest set of connected works is mapped; separate small groups would only stretch the map.
    const adjacency = new Map();
    for (const key of kept.keys()) { const [a, b] = key.split('|'); for (const [u, v] of [[a, b], [b, a]]) { if (!adjacency.has(u)) adjacency.set(u, []); adjacency.get(u).push(v); } }
    let largest = new Set();
    const seen = new Set();
    for (const start of chosen) if (adjacency.has(start) && !seen.has(start)) {
        const component = new Set([start]), queue = [start]; seen.add(start);
        while (queue.length) for (const v of adjacency.get(queue.pop())) if (!seen.has(v)) { seen.add(v); component.add(v); queue.push(v); }
        if (component.size > largest.size) largest = component;
    }
    for (const key of [...kept.keys()]) if (!largest.has(key.split('|')[0])) kept.delete(key);
    const ids = chosen.filter(id => largest.has(id)), index = new Map(ids.map((id, i) => [id, i]));
    const stats = { kind, candidates: candidates.length, linked_works: tls.size, links: links.size, selected: [...adjacency.keys()].length, outside_largest_set: [...adjacency.keys()].length - ids.length };
    if (ids.length < 2) return { kind, nodes: [], edges: [], clusters: [], stats };
    const { strength } = associationStrength(kept);
    const edges = [...kept].map(([key, n]) => { const [a, b] = key.split('|'); return [index.get(a), index.get(b), strength.get(key), n]; });
    const raw = louvain(ids.length, edges.map(([i, j, s]) => [i, j, s]));
    // Clusters are numbered by size as in VOSviewer; clusters smaller than the minimum share cluster 0 ("other").
    const sizes = new Map(); for (const c of raw) sizes.set(c, (sizes.get(c) || 0) + 1);
    const order = [...sizes].filter(([, s]) => s >= MAP_LIMITS.minClusterSize).sort((a, b) => b[1] - a[1] || a[0] - b[0]).map(([c]) => c);
    const clusterOf = raw.map(c => order.indexOf(c) + 1);
    const xy = vosLayout(ids.length, edges.map(([i, j, s]) => [i, j, s]));
    const degree = new Map(), strengthOf = new Map();
    for (const [i, j, , n] of edges) for (const v of [i, j]) { degree.set(v, (degree.get(v) || 0) + 1); strengthOf.set(v, (strengthOf.get(v) || 0) + n); }
    const nodes = ids.map((id, i) => { const w = byId.get(id); return { id, label: label(w), title: w.title, year: w.year || null, citations: w.citations || 0,
        x: Math.round(xy[i][0] * 1e4) / 1e4, y: Math.round(xy[i][1] * 1e4) / 1e4, cluster: clusterOf[i], links: degree.get(i) || 0, strength: strengthOf.get(i) || 0, topic: w.topic || null }; });
    const shown = edges.sort((a, b) => b[2] - a[2] || b[3] - a[3]).slice(0, MAP_LIMITS.edges).map(([i, j, s, n]) => ({ source: ids[i], target: ids[j], strength: n, association: s }));
    const clusters = [...new Set(clusterOf)].sort((a, b) => a - b).map(c => {
        const members = nodes.filter(n => n.cluster === c), topics = new Map();
        for (const n of members) if (n.topic) topics.set(n.topic, (topics.get(n.topic) || 0) + 1);
        const topic = [...topics].sort((a, b) => b[1] - a[1] || (a[0] < b[0] ? -1 : 1))[0]?.[0] || null;
        return { id: c, size: members.length, topic, citations: members.reduce((s, n) => s + n.citations, 0) };
    });
    return { kind, nodes, edges: shown, clusters, stats: { ...stats, nodes: nodes.length, shown_links: shown.length, total_links: edges.length } };
}

// VOSviewer map and network files (tab-separated, numeric ids), readable by VOSviewer's "Create map based on network data".
export function vosviewerFiles(map) {
    const num = new Map(map.nodes.map((n, i) => [n.id, i + 1])), clean = s => String(s ?? '').replace(/[\t\r\n"]+/g, ' ').trim();
    const header = 'id\tlabel\tdescription\turl\tx\ty\tcluster\tweight<Links>\tweight<Total link strength>\tweight<Citations>\tscore<Pub. year>';
    const rows = map.nodes.map(n => [num.get(n.id), clean(n.label), clean(n.title), `https://openalex.org/${n.id}`, n.x, n.y, n.cluster || '', n.links, n.strength, n.citations, n.year ?? ''].join('\t'));
    return { map: [header, ...rows].join('\n') + '\n', network: map.edges.map(e => `${num.get(e.source)}\t${num.get(e.target)}\t${e.strength}`).join('\n') + '\n' };
}
