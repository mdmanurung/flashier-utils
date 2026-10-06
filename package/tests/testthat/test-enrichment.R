test_that("enrichment requires explicit engines and preserves universe losses", {
  d <- association_fixture(K=8)
  pathways <- list(one=c("g1","g2","g3"),empty="unknown")
  expect_error(enrich_factors(d$view,pathways),"explicitly")
  empty <- enrich_factors(d$view,list(none="unknown"),"fgsea")
  expect_equal(nrow(empty$table),0L)
  expect_equal(empty$analysis_metadata$status,"no_overlapping_pathways")
  expect_error(enrich_factors(d$view,pathways,"susie"),"unsupported")
  expect_error(enrich_factors(d$view,pathways,"fgsea",feature_map=c(g1="x",g2="x")),"unique")
  expect_error(.check_engine("flashbridge_missing_engine"),"install it explicitly")
})

test_that("ranked fgsea results agree directly with the selected engine", {
  skip_if_not_installed("fgsea")
  d <- association_fixture(K=8)
  pathways <- list(one=c("g1","g2","g3"),two=c("g4","g5","g6"))
  # Distinct effects avoid native ties in this parity fixture.
  d$view$effects <- matrix(seq_len(64),8,dimnames=list(paste0("g",1:8),paste0("F",1:8)))
  d$view <- .make_representation(d$view$activity,d$view$effects,d$view$manifest)
  settings <- list(seed=4,nPermSimple=100L,eps=1e-3,minSize=2L,scoreType="pos")
  out <- enrich_factors(d$view,pathways,"fgsea",factors="F1",control=settings)
  set.seed(4)
  native <- as.data.frame(fgsea::fgseaMultilevel(pathways,setNames(d$view$effects[,1],rownames(d$view$effects)),BPPARAM=BiocParallel::SerialParam(),nPermSimple=100L,eps=1e-3,minSize=2L,scoreType="pos"))
  expect_equal(out$table$NES,native$NES)
  expect_equal(out$table$pval,native$pval)
  expect_equal(out$analysis_metadata$universe,rownames(d$view$effects))
})

test_that("decoupleR returns native method-specific scores", {
  skip_if_not_installed("decoupleR")
  d <- association_fixture(K=8)
  d$view$effects <- matrix(seq_len(64),8,dimnames=list(paste0("g",1:8),paste0("F",1:8)))
  d$view <- .make_representation(d$view$activity,d$view$effects,d$view$manifest)
  sets <- list(one=c("g1","g2","g3"),two=c("g4","g5","g6"))
  out <- enrich_factors(d$view,sets,"decoupleR",control=list(method="ulm",minsize=2L))
  network <- data.frame(source=rep(names(sets),each=3),target=unlist(sets),mor=1)
  native <- as.data.frame(decoupleR::run_ulm(d$view$effects,network,minsize=2L))
  expect_equal(out$table$score,native$score)
  expect_equal(out$analysis_metadata$interpretation,"network_activity_score")
  mlm <- enrich_factors(d$view,sets,"decoupleR",control=list(method="mlm",minsize=2L))
  direct <- as.data.frame(decoupleR::run_mlm(d$view$effects,network,minsize=2L))
  expect_equal(mlm$table$score,direct$score)
  expect_equal(mlm$table$p_value,direct$p_value)
})

test_that("weighted network annotations preserve native scores and mapped universes", {
  skip_if_not_installed("decoupleR")
  set.seed(31)
  d <- association_fixture(K=8)
  d$view <- .make_representation(d$view$activity,
    matrix(rnorm(64),8,dimnames=dimnames(d$view$effects)),d$view$manifest)
  network <- data.frame(source=rep(c("TF1","pathway1"),each=4),
    target=c("a","b","c","absent","d","e","f","absent"),
    mor=c(1,-2,0.5,1,-1,2,0.3,-2))
  attr(network,"organism") <- "human"
  attr(network,"annotation_version") <- "test_v1"
  mapping <- setNames(letters[1:7],paste0("g",1:7))
  B <- d$view$effects[1:7,c("F2","F1")]
  rownames(B) <- letters[1:7]
  for (method in c("ulm","mlm")) {
    out <- enrich_factors(d$view,network,"decoupleR",factors=c("F2","F1"),
      feature_map=mapping,control=list(method=method,minsize=2L))
    native <- as.data.frame(getExportedValue("decoupleR",paste0("run_",method))(
      B,network[network$target != "absent",],minsize=2L))
    native <- native[order(match(native$condition,colnames(B))),,drop=FALSE]
    rownames(native) <- NULL
    expect_equal(out$table[,names(native)],native,ignore_attr=TRUE)
    expect_equal(out$analysis_metadata$universe,letters[1:7])
    expect_equal(out$analysis_metadata$mapping_excluded,"g8")
    expect_equal(out$analysis_metadata$network_targets_excluded,"absent")
    expect_equal(out$analysis_metadata$pathway_sizes$overlap_size,c(3L,3L))
    expect_equal(out$analysis_metadata$organism,"human")
    expect_equal(out$analysis_metadata$annotation_version,"test_v1")
    expect_true(all(out$table$basis_id == d$view$manifest$basis_id))
    expect_equal(unique(out$table$condition), c("F2", "F1"))
    single <- enrich_factors(d$view,network,"decoupleR",factors="F2",
      feature_map=mapping,control=list(method=method,minsize=2L))
    expect_equal(single$table$score,out$table$score[out$table$condition == "F2"])
  }
  flipped <- .make_representation(-d$view$activity,-d$view$effects,d$view$manifest)
  original <- enrich_factors(d$view,network,"decoupleR",feature_map=mapping,control=list(minsize=2L))
  reverse <- enrich_factors(flipped,network,"decoupleR",feature_map=mapping,control=list(minsize=2L))
  expect_equal(reverse$table$score,-original$table$score)
  expect_equal(reverse$table$p_value,original$table$p_value)
  negative <- network; negative$mor <- -negative$mor
  reverse_weights <- enrich_factors(d$view,negative,"decoupleR",feature_map=mapping,control=list(minsize=2L))
  expect_equal(reverse_weights$table$score,-original$table$score)
  expect_equal(reverse_weights$table$p_value,original$table$p_value)
  expect_error(decoupleR::run_ulm(B,network,minsize=9L))
  expect_error(enrich_factors(d$view,network,"decoupleR",feature_map=mapping,control=list(minsize=9L)))
  network$mor[1] <- -network$mor[1]
  changed <- enrich_factors(d$view,network,"decoupleR",feature_map=mapping,control=list(minsize=2L))
  expect_false(identical(changed$analysis_metadata$network_hash,original$analysis_metadata$network_hash))
})

test_that("network validation fails before native dispatch", {
  view <- fixture_view()
  net <- data.frame(source=c("TF1","TF1"),target=c("g1","g2"),mor=c(1,-1))
  expect_error(enrich_factors(view,net,"fgsea"),"Weighted networks require")
  expect_error(enrich_factors(view,net[,1:2],"decoupleR"),"mor columns")
  bad <- net; bad$mor[1] <- Inf
  expect_error(enrich_factors(view,bad,"decoupleR"),"finite numeric")
  bad <- net; bad$source[1] <- NA_character_
  expect_error(enrich_factors(view,bad,"decoupleR"),"nonempty strings")
  bad <- net; bad$target[1] <- " "
  expect_error(enrich_factors(view,bad,"decoupleR"),"nonempty strings")
  bad <- net; bad$mor <- c("1","-1")
  expect_error(enrich_factors(view,bad,"decoupleR"),"finite numeric")
  bad <- net; bad$mor[1] <- NaN
  expect_error(enrich_factors(view,bad,"decoupleR"),"finite numeric")
  expect_error(enrich_factors(view,rbind(net,net[1,]),"decoupleR"),"edges must be unique")
  net$target <- c("absent1","absent2")
  out <- enrich_factors(view,net,"decoupleR")
  expect_equal(nrow(out$table),0L)
  expect_equal(out$analysis_metadata$network_targets_excluded,net$target)
})
