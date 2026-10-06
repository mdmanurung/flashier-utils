test_that("prediction R squared is negative even with perfect correlation", {
  observed <- data.frame(feature=c("g1","g2","g3"),estimate=c(1,2,3),contrast="A - B")
  reconstructed <- observed; reconstructed$estimate <- c(11,12,13)
  out <- validate_backprojection(reconstructed,observed,feature_scale="arcsinh_intensity",top_n=20L,plot=FALSE)
  expect_equal(out$metrics$value[out$metrics$metric=="prediction_r2"],-149)
  expect_equal(out$metrics$value[out$metrics$metric=="correlation_squared"],1)
  expect_equal(out$top_n_overlap$effective_n,3)
  expect_equal(out$residuals$residual,c(-10,-10,-10))
  attr(reconstructed,"analysis_metadata") <- list(feature_scale="log1p")
  expect_error(validate_backprojection(reconstructed,observed,feature_scale="arcsinh_intensity"),"scales")
  attr(reconstructed,"analysis_metadata") <- NULL; reconstructed$contrast <- "B - A"
  expect_error(validate_backprojection(reconstructed,observed,feature_scale="arcsinh_intensity"),"Contrast")
  observed$estimate <- 1
  out <- validate_backprojection(observed,observed,feature_scale="arcsinh_intensity",plot=FALSE)
  expect_true(is.na(out$metrics$value[out$metrics$metric=="prediction_r2"]))
})
