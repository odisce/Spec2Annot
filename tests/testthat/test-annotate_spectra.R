test_that(
  "Spec2Annot::annotat_spectra", {
    spectra_in <- Spec2Annot::spectra_full
    ## Create diff grid
    # Rcpp::sourceCpp("./src/cpp_annotate_spectra.cpp")
    annotate_spectra(spectra_in[sample(seq_len(.N), 10), mz])
    find_comb(10, 2)
  }
)