context("test-produce_heatmap")

valid_cov <- import_bedgraphs(get_demo_bdg())
valid_narrow_peak_file <- peakimport::get_demo_file("bed")
valid_peaks <- peakimport::import_peaks(valid_narrow_peak_file)
valid_peaks_chr19 <- valid_peaks[seqnames(valid_peaks) == "chr19"]

test_that("cov works as expected", {
    # Valid case
    heatmap <- produce_heatmap(cov = valid_cov,
                               peaks = valid_peaks_chr19,
                               name = "test")
    expect_is(heatmap, "Heatmap")
    
    # Invalid cases
    expect_error(produce_heatmap(cov = GRanges(),
                                 peaks = valid_peaks_chr19,
                                 name = "test"),
                 "length(cov) > 0 is not TRUE",
                 fixed = TRUE)
    expect_error(produce_heatmap(cov = "",
                                 peaks = valid_peaks_chr19,
                                 name = "test"),
                 "is(cov, \"GRanges\") is not TRUE",
                 fixed = TRUE)
    invalid_cov_style <- valid_cov
    seqlevelsStyle(invalid_cov_style) <- "NCBI"
    expect_error(produce_heatmap(cov = invalid_cov_style,
                                 peaks = valid_peaks_chr19,
                                 name = "test"),
                 "all(seqnames_peaks %in% seqnames_cov) is not TRUE",
                 fixed = TRUE)
})

test_that("peaks works as expected", {
    # Valid case
    ## Tested in "cov works as expected tests"

    # Invalid cases
    expect_error(produce_heatmap(cov = valid_cov,
                                 peaks = GRanges(),
                                 name = "test"),
                 "length(peaks) > 0 is not TRUE",
                 fixed = TRUE)
    expect_error(produce_heatmap(cov = valid_cov,
                                 peaks = "",
                                 name = "test"),
                 "is(peaks, \"GRanges\") is not TRUE",
                 fixed = TRUE)
    expect_error(produce_heatmap(cov = valid_cov,
                                 peaks = valid_peaks,
                                 name = "test"),
                 "all(seqnames_peaks %in% seqnames_cov) is not TRUE",
                 fixed = TRUE)
})

test_that("force_seqlevels works as expected", {
    # Valid case
    heatmap <- produce_heatmap(cov = valid_cov,
                               peaks = valid_peaks,
                               name = "test",
                               force_seqlevels = TRUE)
    expect_is(heatmap, "Heatmap")

    # Invalid cases
    expect_error(produce_heatmap(cov = valid_cov,
                                 peaks = valid_peaks,
                                 name = "test",
                                 force_seqlevels = FALSE),
                 "all(seqnames_peaks %in% seqnames_cov) is not TRUE",
                 fixed = TRUE)
    expect_error(produce_heatmap(cov = valid_cov,
                                 peaks = valid_peaks_chr19,
                                 name = "test",
                                 force_seqlevels = ""),
                 "is(force_seqlevels, \"logical\") is not TRUE",
                 fixed = TRUE)
})
