

collapse_mask <- function(helpers = TRUE) {

  detach("package:collapse", character.only = TRUE)
  # `%=%`(c("subset", "transform"), list(fsubset, ftransform), env = ns)
  ns <- parent.env(environment()) # getNamespace('collapse')
  # ns$subset <- ns$fsubset
  # assign("subset", ns$fsubset, envir = ns)
  subset <- ns$fsubset
  environment(subset) <- ns
  namespaceExport(ns, "subset")
  # assignInMyNamespace("subset", ns$fsubset)

  # fixInNamespace("transform", "stats")

  attachNamespace(ns)

  # library(collapse)

  # assignInMyNamespace("subset", fsubset)
  # assignInMyNamespace("transform", ftransform)
  # assignInMyNamespace("transform<-", `ftransform<-`)
  # assignInMyNamespace("transformv", ftransformv)
  # assignInMyNamespace("compute", fcompute)
  # assignInMyNamespace("computev", fcomputev)
  # assignInMyNamespace("select", fselect)
  # assignInMyNamespace("select<-", `fselect<-`)
  # assignInMyNamespace("group_by", fgroup_by)
  # assignInMyNamespace("group_vars", fgroup_vars)
  # assignInMyNamespace("ungroup", fungroup)
  # assignInMyNamespace("summarise", fsummarise)
  # assignInMyNamespace("mutate", fmutate)
  # assignInMyNamespace("rename", frename)
  #
  # if(helpers) {
  #   assignInMyNamespace("droplevels", fdroplevels)
  #   assignInMyNamespace("interaction", finteraction)
  #   assignInMyNamespace("nlevels", fnlevels)
  #   assignInMyNamespace("unique", funique)
  #   assignInMyNamespace("nrow", fnrow)
  #   assignInMyNamespace("ncol", fncol)
  #   assignInMyNamespace("dim", fdim)
  # }
}

collapse_F_to_FALSE <- function() assignInMyNamespace("F", FALSE)
