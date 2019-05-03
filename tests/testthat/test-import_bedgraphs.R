context("test-import_bedgraphs")

valid_filename <- get_demo_bdg()
valid_filename_idx <- get_demo_bdg_idx()

all_seqlevels_plus_invalid <-
    system.file("extdata/invalid_all_seqlevels_plus_invalid.bdg",
                package = "heatmaps", mustWork = TRUE)
invalid_bdg_format <- system.file("extdata/invalid_bdg_format.txt",
                                  package = "heatmaps",
                                  mustWork = TRUE)
one_seqlevel_plus_invalid <-
    system.file("extdata/one_seqlevel_plus_invalid.bdg",
                package = "heatmaps", mustWork = TRUE)

valid_gr <- function(gr, expected_length = 30000) {
    expect_is(gr, "GRanges")
    expect_equal(length(gr), expected_length)
}

test_that("filename works as expected", {
    # Valid cases
    valid_gr(import_bedgraphs(filename = valid_filename))
    valid_gr(import_bedgraphs(filename = valid_filename))
    
    # Invalid cases
    expect_error(import_bedgraphs(filename = "filename_not_found"),
                 "file.exists(filename) is not TRUE",
                 fixed = TRUE)
    expect_error(import_bedgraphs(filename = invalid_bdg_format))
    expect_error(import_bedgraphs(filename = TRUE))
})

test_that("filter_negative_coverages works as expected", {
    # Valid cases
    valid_gr(import_bedgraphs(filename = valid_filename,
                              filter_negative_coverages = TRUE))
    valid_gr(import_bedgraphs(filename = valid_filename,
                              filter_negative_coverages = FALSE))

    # Invalid case
    expect_error(import_bedgraphs(filename = valid_filename,
                                  filter_negative_coverages = "TRUE"),
                 "is.logical(filter_negative_coverages) is not TRUE",
                 fixed = TRUE)
})

test_that("keep_standard_chromosomes works as expected", {
    # Valid cases
    valid_gr(import_bedgraphs(filename = valid_filename,
                              keep_standard_chromosomes = TRUE))
    gr <- import_bedgraphs(filename = valid_filename,
                           keep_standard_chromosomes = FALSE)
    expect_is(gr, "GRanges")
    expect_equal(length(gr), 30000+2)

    # Invalid case
    expect_error(import_bedgraphs(filename = valid_filename,
                                  keep_standard_chromosomes = "TRUE"),
                 "is.logical(keep_standard_chromosomes) is not TRUE",
                 fixed = TRUE)
})

test_that("genome works as expected", {
    # Valid case
    gr <- import_bedgraphs(filename = valid_filename, genome = "hg38")
    valid_gr(gr)
    expect_equal(unique(as.data.frame(seqinfo(gr))$genome), "hg38")

    gr <- import_bedgraphs(filename = all_seqlevels_plus_invalid,
                           genome = "hg38")
    valid_gr(gr, 25)
    expect_equal(unique(as.data.frame(seqinfo(gr))$genome), "hg38")

    # Invalid case
    expect_error(import_bedgraphs(filename = valid_filename,
                                  genome = "invalid_genome"),
                 "any(names(si) == genome) is not TRUE",
                 fixed = TRUE)

    expect_error(import_bedgraphs(filename = all_seqlevels_plus_invalid,
                                 genome = "hg38",
								 keep_standard_chromosomes = FALSE))
    expect_error(import_bedgraphs(filename = one_seqlevel_plus_invalid,
                                 genome = "hg38",
								 keep_standard_chromosomes = FALSE))
})
