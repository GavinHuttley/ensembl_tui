BUG

- Installing an alignment no longer writes duplicate rows. Ensembl emits an
  alignment block once for each segment of the reference species it holds, so
  the same record reaches `add_records()` several times in one batch. The
  existing check only skipped block ids written by an earlier call, so those
  copies were all stored: 5,490 of 60,017 rows (9.1%) in a three species
  primate install. Query results were already correct, because the read path
  collapses identical records into a set, so an existing installed store does
  not need rebuilding.

- `AlignRecord.__eq__()` raised `ValueError` when comparing two records with
  the same coordinates but gap arrays of different shapes. Such records hash
  alike, so they met in the set built by `AlignDb.get_records_matching()`.
  They now compare with `numpy.array_equal()`, which tolerates a shape
  mismatch.
