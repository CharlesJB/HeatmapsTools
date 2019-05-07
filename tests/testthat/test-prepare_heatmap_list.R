context("test-prepare_heatmap_list")

htl_w_partition <- get_demo_heatmap_list(partitions = TRUE)
htl_wo_partition <- get_demo_heatmap_list(partitions = FALSE)

test_that("heatmap_list works as expected", {
    # Valid case
    heatmap_list <- prepare_heatmap_list(heatmap_list = htl_wo_partition)
    expect_is(heatmap_list, "HeatmapList")

    # Invalid cases
    expect_error(prepare_heatmap_list(heatmap_list = ""),
                 "is(heatmap_list, \"list\") is not TRUE",
                 fixed = TRUE)
    expect_error(prepare_heatmap_list(heatmap_list = ""),
                 "is(heatmap_list, \"list\") is not TRUE",
                 fixed = TRUE)
    expect_error(prepare_heatmap_list(heatmap_list = list(a = "a")),
                 "all(map_lgl(heatmap_list, is, \"Heatmap\")) is not TRUE",
                 fixed = TRUE)
    expect_error(prepare_heatmap_list(heatmap_list = c(htl_wo_partition, a = "a")),
                 "all(map_lgl(heatmap_list, is, \"Heatmap\")) is not TRUE",
                 fixed = TRUE)
})

test_that("partitions works as expected", {
    partitions <- as.character(c(rep(1,10), rep(2,10)))

    # valid cases
    heatmap_list <- prepare_heatmap_list(heatmap_list = htl_w_partition,
                                         partitions = partitions)
    expect_is(heatmap_list, "HeatmapList")

    # Invalid cases
    expect_error(prepare_heatmap_list(heatmap_list = htl_w_partition,
                                      partitions = as.numeric(partitions)),
                 "is.character(partitions) is not TRUE",
                 fixed = TRUE)
    expect_error(prepare_heatmap_list(heatmap_list = htl_w_partition,
                                      partitions = ""),
                 "length(partitions) == nrow(heatmap_list[[1]]@matrix) is not TRUE",
                 fixed = TRUE)
    expect_error(prepare_heatmap_list(heatmap_list = htl_w_partition,
                                      partitions = rep("1", 20)),
                 "length(unique(partitions)) > 1 is not TRUE",
                 fixed = TRUE)
})
