use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn get_motif_with_idx(indices: &[usize], peptide: &str) -> Vec<u8> {
    let pep_len = peptide.len();
    indices
        .iter()
        .flat_map(|idx| {
            if *idx < pep_len {
                Some(peptide.as_bytes()[*idx])
            } else {
                None
            }
        })
        .collect()
}

fn get_motif_with_iter(indices: &[usize], peptide: &str) -> Vec<u8> {
    indices
        .iter()
        .flat_map(|idx| peptide.bytes().nth(*idx))
        .collect()
}
fn criterion_benchmark(c: &mut Criterion) {
    let input_peptide = black_box("ABCDEFGHIJ");
    let indices = black_box(vec![2, 3, 4, 5, 6, 9]);

    c.bench_function("motif_with_idx 9mer CD8", |b| {
        b.iter(|| get_motif_with_idx(&indices[..], &input_peptide))
    });

    c.bench_function("motif_with_iter 9mer CD8", |b| {
        b.iter(|| get_motif_with_iter(&indices[..], &input_peptide))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
