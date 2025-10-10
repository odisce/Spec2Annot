test_that(
  "Spec2Annot::annotat_spectra", {
    spectra_in <- Spec2Annot::spectra_full
    ## Create diff grid
    # Rcpp::sourceCpp("./src/cpp_annotate_spectra.cpp")
    res <- annotate_spectra(spectra_in[sample(seq_len(.N), 10), mz] %>% sort())
    find_comb(10, 2)
    
    # devtools::load_all()
    # combn_cpp(spectra_in[, unique(mz)][1:10], 2)
  }
)