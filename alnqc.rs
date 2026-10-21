// Rust implementations of alnqc.

#[allow(non_snake_case)]
#[derive(Clone, Debug, Default)]
pub struct AqRefs {
    pub lengths: Vec<i64>,
    pub genes: Vec<i64>,
    pub geneIds: Vec<String>,
    pub biotypes: Vec<i64>,
    pub biotypeNames: Vec<String>,
}

#[derive(Clone, Debug, Default)]
pub struct AqTxPoints {
    pub tx: i64,
    pub fragments: i64,
    pub points: Vec<i64>,
}

#[allow(non_snake_case)]
#[derive(Clone, Debug, Default)]
pub struct AqAlnQc {
    pub fragments: i64,
    pub unmapped: i64,
    pub decoy: i64,
    pub oneTx: i64,
    pub oneGene: i64,
    pub multiGene: i64,
    pub biotypeCounts: Vec<f64>,
    pub txPoints: Vec<AqTxPoints>,
}

const AQ_POINTS: usize = 100;
const AQ_MIN_TX_LEN: i64 = 100;
const AQ_TIN_MIN_FRAGMENTS: i64 = 10;

pub fn aq_refs(names: &Vec<String>, lengths: &Vec<i64>, txs: &Vec<BioTxGene>) -> AqRefs {
    if names.len() != lengths.len() {
        panic!("alnQcRefs: {} names but {} lengths", names.len(), lengths.len());
    }
    let by_tx: std::collections::HashMap<&str, &BioTxGene> = txs.iter().map(|t| (t.tx.as_str(), t)).collect();
    let mut r = AqRefs { lengths: lengths.clone(), ..Default::default() };
    let mut gene_no: std::collections::HashMap<String, i64> = std::collections::HashMap::new();
    let mut type_no: std::collections::HashMap<String, i64> = std::collections::HashMap::new();
    for n in names {
        match by_tx.get(n.as_str()) {
            None => {
                r.genes.push(-1);
                r.biotypes.push(-1);
            }
            Some(t) => {
                let g = *gene_no.entry(t.gene.clone()).or_insert_with(|| {
                    r.geneIds.push(t.gene.clone());
                    r.geneIds.len() as i64 - 1
                });
                let bt = if t.biotype.is_empty() { "unknown".to_string() } else { t.biotype.clone() };
                let b = *type_no.entry(bt.clone()).or_insert_with(|| {
                    r.biotypeNames.push(bt);
                    r.biotypeNames.len() as i64 - 1
                });
                r.genes.push(g);
                r.biotypes.push(b);
            }
        }
    }
    r
}

// The sample position of point k (0-based) on a transcript of length `len`.
fn aq_point(len: i64, k: usize) -> i64 {
    (len - 1) * (k as i64 + 1) / AQ_POINTS as i64
}

fn aq_cover(points: &mut [i64], len: i64, a: i64, b: i64) {
    let (a, b) = (a.max(0), b.min(len));
    if a >= b {
        return;
    }
    for k in 0..AQ_POINTS {
        if aq_point(len, k) < a {
            continue;
        }
        if aq_point(len, k) >= b {
            break;
        }
        points[k] += 1;
    }
}

// Mate lengths: (read 1, read 2); read 2 is 0 for single-end reads.
fn aq_chunk(refs: &AqRefs, lens: &[(i64, i64)], xs: &[BioFragHits]) -> AqAlnQc {
    if lens.len() != xs.len() {
        panic!("alnQc: {} reads but {} placements", lens.len(), xs.len());
    }
    let mut q = AqAlnQc { biotypeCounts: vec![0.0; refs.biotypeNames.len()], ..Default::default() };
    let mut pts: std::collections::HashMap<i64, AqTxPoints> = std::collections::HashMap::new();
    for (f, &(l1, l2)) in xs.iter().zip(lens) {
        q.fragments += 1;
        if f.hits.is_empty() {
            if f.decoy { q.decoy += 1 } else { q.unmapped += 1 }
            continue;
        }
        let gene_of = |t: i64| refs.genes.get(t as usize).copied().unwrap_or(-1);
        let tx0 = f.hits[0].tx;
        let g0 = gene_of(tx0);
        let one_tx = f.hits.iter().all(|h| h.tx == tx0);
        let one_gene = one_tx || (g0 >= 0 && f.hits.iter().all(|h| gene_of(h.tx) == g0));
        if one_tx {
            q.oneTx += 1;
        } else if one_gene {
            q.oneGene += 1;
        } else {
            q.multiGene += 1;
        }
        let total: f64 = f.hits.iter().map(|h| h.weight.max(0.0)).sum();
        for h in &f.hits {
            let b = refs.biotypes.get(h.tx as usize).copied().unwrap_or(-1);
            if b >= 0 {
                let w = if total > 0.0 { h.weight.max(0.0) / total } else { 1.0 / f.hits.len() as f64 };
                q.biotypeCounts[b as usize] += w;
            }
        }
        if !one_gene {
            continue;
        }
        let mut seen: Vec<i64> = Vec::new();
        for h in &f.hits {
            let len = refs.lengths.get(h.tx as usize).copied().unwrap_or(0);
            if len < AQ_MIN_TX_LEN {
                continue;
            }
            let e = pts.entry(h.tx).or_insert_with(|| AqTxPoints { tx: h.tx, fragments: 0, points: vec![0; AQ_POINTS] });
            if !seen.contains(&h.tx) {
                e.fragments += 1;
                seen.push(h.tx);
            }
            match h.mates {
                0 if h.fragLen > 0 => {
                    let (left, right) = if h.fwd { (l1, l2) } else { (l2, l1) };
                    aq_cover(&mut e.points, len, h.pos, h.pos + left);
                    let rs = (h.pos + h.fragLen - right).max(h.pos + left);
                    aq_cover(&mut e.points, len, rs, h.pos + h.fragLen);
                }
                2 => aq_cover(&mut e.points, len, h.pos, h.pos + l2),
                _ => aq_cover(&mut e.points, len, h.pos, h.pos + l1),
            }
        }
    }
    let mut v: Vec<AqTxPoints> = pts.into_values().collect();
    v.sort_by_key(|p| p.tx);
    q.txPoints = v;
    q
}

pub fn aq_aln_qc_pairs(_anchor: i64, refs: &AqRefs, reads: &Vec<(BioFastqEntry, BioFastqEntry)>, hits: &Vec<BioFragHits>) -> AqAlnQc {
    let lens: Vec<(i64, i64)> = reads.iter().map(|(a, b)| (a.bases.bases.len() as i64, b.bases.bases.len() as i64)).collect();
    aq_chunk(refs, &lens, hits)
}

pub fn aq_aln_qc_reads(_anchor: i64, refs: &AqRefs, reads: &Vec<BioFastqEntry>, hits: &Vec<BioFragHits>) -> AqAlnQc {
    let lens: Vec<(i64, i64)> = reads.iter().map(|a| (a.bases.bases.len() as i64, 0)).collect();
    aq_chunk(refs, &lens, hits)
}

pub fn aq_merge(a: &AqAlnQc, b: &AqAlnQc) -> AqAlnQc {
    let n = a.biotypeCounts.len().max(b.biotypeCounts.len());
    let mut bt = vec![0.0; n];
    for (i, x) in a.biotypeCounts.iter().enumerate() {
        bt[i] += x;
    }
    for (i, x) in b.biotypeCounts.iter().enumerate() {
        bt[i] += x;
    }
    let mut pts = Vec::with_capacity(a.txPoints.len() + b.txPoints.len());
    let (mut i, mut j) = (0, 0);
    while i < a.txPoints.len() || j < b.txPoints.len() {
        let take_a = j >= b.txPoints.len() || (i < a.txPoints.len() && a.txPoints[i].tx < b.txPoints[j].tx);
        let take_b = i >= a.txPoints.len() || (j < b.txPoints.len() && b.txPoints[j].tx < a.txPoints[i].tx);
        if take_a {
            pts.push(a.txPoints[i].clone());
            i += 1;
        } else if take_b {
            pts.push(b.txPoints[j].clone());
            j += 1;
        } else {
            let (x, y) = (&a.txPoints[i], &b.txPoints[j]);
            pts.push(AqTxPoints { tx: x.tx, fragments: x.fragments + y.fragments, points: x.points.iter().zip(&y.points).map(|(p, q)| p + q).collect() });
            i += 1;
            j += 1;
        }
    }
    AqAlnQc {
        fragments: a.fragments + b.fragments,
        unmapped: a.unmapped + b.unmapped,
        decoy: a.decoy + b.decoy,
        oneTx: a.oneTx + b.oneTx,
        oneGene: a.oneGene + b.oneGene,
        multiGene: a.multiGene + b.multiGene,
        biotypeCounts: bt,
        txPoints: pts,
    }
}

pub fn aq_gene_body_coverage(q: &AqAlnQc) -> Vec<i64> {
    let mut out = vec![0i64; AQ_POINTS];
    for p in &q.txPoints {
        for (o, x) in out.iter_mut().zip(&p.points) {
            *o += x;
        }
    }
    out
}

fn aq_tin(points: &[i64]) -> f64 {
    let total: i64 = points.iter().sum();
    if total == 0 {
        return 0.0;
    }
    let h: f64 = points.iter().filter(|&&c| c > 0).map(|&c| {
        let p = c as f64 / total as f64;
        -p * p.ln()
    }).sum();
    100.0 * h.exp() / AQ_POINTS as f64
}

pub fn aq_tin_values(min_fragments: i64, q: &AqAlnQc) -> Vec<(i64, f64)> {
    q.txPoints.iter().filter(|p| p.fragments >= min_fragments).map(|p| (p.tx, aq_tin(&p.points))).collect()
}

pub fn aq_gene_body_coverage_text(sample: &String, q: &AqAlnQc) -> String {
    let cov = aq_gene_body_coverage(q);
    let head: Vec<String> = (1..=AQ_POINTS).map(|i| i.to_string()).collect();
    let vals: Vec<String> = cov.iter().map(|c| format!("{}.0", c)).collect();
    format!("Percentile\t{}\n{}\t{}\n", head.join("\t"), sample, vals.join("\t"))
}

pub fn aq_tin_summary_text(sample: &String, q: &AqAlnQc) -> String {
    let mut t: Vec<f64> = aq_tin_values(AQ_TIN_MIN_FRAGMENTS, q).into_iter().map(|x| x.1).collect();
    let (mean, median, sd) = if t.is_empty() {
        (0.0, 0.0, 0.0)
    } else {
        t.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let n = t.len() as f64;
        let mean = t.iter().sum::<f64>() / n;
        let median = if t.len() % 2 == 1 { t[t.len() / 2] } else { (t[t.len() / 2 - 1] + t[t.len() / 2]) / 2.0 };
        let sd = (t.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n).sqrt();
        (mean, median, sd)
    };
    format!("Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)\n{}\t{}\t{}\t{}\n", sample, mean, median, sd)
}

pub fn aq_biotype_counts_text(refs: &AqRefs, q: &AqAlnQc) -> String {
    let mut s = String::from("# id: 'biotype_counts'\n# section_name: 'Biotype Counts'\n# description: \"shows fragments placed on transcripts of each gene biotype, split by placement weight.\"\n# plot_type: 'bargraph'\n# anchor: 'featurecounts_biotype'\n# pconfig:\n#     id: \"featurecounts_biotype_plot\"\n#     title: \"Biotypes\"\n#     xlab: \"# Reads\"\n#     cpswitch_counts_label: \"Number of Reads\"\n");
    for (n, c) in refs.biotypeNames.iter().zip(&q.biotypeCounts) {
        s.push_str(&format!("{}\t{}\n", n, c));
    }
    s
}

pub fn aq_rrna_percent_text(sample: &String, refs: &AqRefs, q: &AqAlnQc) -> String {
    let total: f64 = q.biotypeCounts.iter().sum();
    let rrna: f64 = refs.biotypeNames.iter().zip(&q.biotypeCounts).filter(|(n, _)| n.as_str() == "rRNA").map(|x| *x.1).sum();
    let pct = if total > 0.0 { 100.0 * rrna / total } else { 0.0 };
    format!("#id: 'biotype-gs'\n#plot_type: 'generalstats'\n#pconfig:\n#    percent_rRNA:\n#        title: '% rRNA'\n#        namespace: 'Biotype Counts'\n#        description: '% reads overlapping rRNA features'\n#        max: 100\n#        min: 0\n#        scale: 'RdYlGn-rev'\nSample\tpercent_rRNA\n'{}'\t{}\n", sample, pct)
}

pub fn aq_read_origin_text(sample: &String, q: &AqAlnQc) -> String {
    let n = q.fragments.max(1) as f64;
    let pct = |x: i64| 100.0 * x as f64 / n;
    format!(
        "#id: 'read_origin'\n#section_name: 'Read origin'\n#description: 'Fragments by where the mapper placed them.'\n#plot_type: 'bargraph'\n#pconfig:\n#    id: 'read_origin_plot'\n#    title: 'Read origin'\n#    cpswitch_c_active: False\nSample\tOne transcript\tOne gene\tSeveral genes\tDecoy (genomic)\tUnmapped\n{}\t{}\t{}\t{}\t{}\t{}\n",
        sample, pct(q.oneTx), pct(q.oneGene), pct(q.multiGene), pct(q.decoy), pct(q.unmapped)
    )
}
