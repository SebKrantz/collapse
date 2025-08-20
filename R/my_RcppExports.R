
BWCpp <- function(x, ng = 0L, g = 0L, gs = NULL, w = NULL, narm = TRUE, theta = 1, set_mean = 0, B = FALSE, fill = FALSE) {
    .Call(Cpp_BW, x, ng, g, gs, w, narm, theta, set_mean, B, fill)
}

BWmCpp <- function(x, ng = 0L, g = 0L, gs = NULL, w = NULL, narm = TRUE, theta = 1, set_mean = 0, B = FALSE, fill = FALSE) {
    .Call(Cpp_BWm, x, ng, g, gs, w, narm, theta, set_mean, B, fill)
}

BWlCpp <- function(x, ng = 0L, g = 0L, gs = NULL, w = NULL, narm = TRUE, theta = 1, set_mean = 0, B = FALSE, fill = FALSE) {
    .Call(Cpp_BWl, x, ng, g, gs, w, narm, theta, set_mean, B, fill)
}

TRAC <- function(x, xAG, g = 0L, ret = 1L, set = FALSE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(set) return(invisible(.Call(C_TRA, x, xAG, g, ret, set)))
  .Call(C_TRA, x, xAG, g, ret, set)
}

TRAmC <- function(x, xAG, g = 0L, ret = 1L, set = FALSE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(set) return(invisible(.Call(C_TRAm, x, xAG, g, ret, set)))
  .Call(C_TRAm, x, xAG, g, ret, set)
}

TRAlC <- function(x, xAG, g = 0L, ret = 1L, set = FALSE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(set) return(invisible(.Call(C_TRAl, x, xAG, g, ret, set)))
  .Call(C_TRAl, x, xAG, g, ret, set)
}

fndistinctC <- function(x, g = NULL, narm = TRUE, nthreads = 1L) {
    .Call(C_fndistinct, x, g, narm, nthreads)
}

pwnobsmCpp <- function(x) {
    .Call(Cpp_pwnobsm, x)
}

fnobsC <- function(x, ng = 0L, g = 0L) {
    .Call(C_fnobs, x, ng, g)
}

varyingCpp <- function(x, ng = 0L, g = 0L, any_group = TRUE) {
    .Call(Cpp_varying, x, ng, g, any_group)
}

varyingmCpp <- function(x, ng = 0L, g = 0L, any_group = TRUE, drop = TRUE) {
    .Call(Cpp_varyingm, x, ng, g, any_group, drop)
}

varyinglCpp <- function(x, ng = 0L, g = 0L, any_group = TRUE, drop = TRUE) {
    .Call(Cpp_varyingl, x, ng, g, any_group, drop)
}

fbstatsCpp <- function(x, ext = FALSE, ng = 0L, g = 0L, npg = 0L, pg = 0L, w = NULL, stable.algo = TRUE, array = TRUE, setn = TRUE, gn = NULL) {
  .Call(Cpp_fbstats, x, ext, ng, g, npg, pg, w, stable.algo, array, setn, gn)
}

fbstatsmCpp <- function(x, ext = FALSE, ng = 0L, g = 0L, npg = 0L, pg = 0L, w = NULL, stable.algo = TRUE, array = TRUE, gn = NULL) {
  .Call(Cpp_fbstatsm, x, ext, ng, g, npg, pg, w, stable.algo, array, gn)
}

fbstatslCpp <- function(x, ext = FALSE, ng = 0L, g = 0L, npg = 0L, pg = 0L, w = NULL, stable.algo = TRUE, array = TRUE, gn = NULL) {
  .Call(Cpp_fbstatsl, x, ext, ng, g, npg, pg, w, stable.algo, array, gn)
}

fdiffgrowthCpp <- function(x, n = 1L, diff = 1L, fill = NA_real_, ng = 0L, g = 0L, gs = NULL, t = NULL, ret = 1L, rho = 1, names = TRUE, power = 1) {
    .Call(Cpp_fdiffgrowth, x, n, diff, fill, ng, g, gs, t, ret, rho, names, power)
}

fdiffgrowthmCpp <- function(x, n = 1L, diff = 1L, fill = NA_real_, ng = 0L, g = 0L, gs = NULL, t = NULL, ret = 1L, rho = 1, names = TRUE, power = 1) {
    .Call(Cpp_fdiffgrowthm, x, n, diff, fill, ng, g, gs, t, ret, rho, names, power)
}

fdiffgrowthlCpp <- function(x, n = 1L, diff = 1L, fill = NA_real_, ng = 0L, g = 0L, gs = NULL, t = NULL, ret = 1L, rho = 1, names = TRUE, power = 1) {
    .Call(Cpp_fdiffgrowthl, x, n, diff, fill, ng, g, gs, t, ret, rho, names, power)
}

flagleadCpp <- function(x, n = 1L, fill = NULL, ng = 0L, g = 0L, t = NULL, names = TRUE) {
    .Call(Cpp_flaglead, x, n, fill, ng, g, t, names)
}

flagleadmCpp <- function(x, n = 1L, fill = NULL, ng = 0L, g = 0L, t = NULL, names = TRUE) {
    .Call(Cpp_flagleadm, x, n, fill, ng, g, t, names)
}

flagleadlCpp <- function(x, n = 1L, fill = NULL, ng = 0L, g = 0L, t = NULL, names = TRUE) {
    .Call(Cpp_flagleadl, x, n, fill, ng, g, t, names)
}


# fnthC <- function(x, n = 0.5, g = NULL, w = NULL, narm = TRUE, ret = 1L, nthreads = 1L, o = NULL, check.o = FALSE) {
#   .Call(C_fnth, x, n, g, w, narm, ret, nthreads, o, check.o)
# }
#
# fnthmC <- function(x, n = 0.5, g = NULL, w = NULL, narm = TRUE, drop = TRUE, ret = 1L, nthreads = 1L) {
#   .Call(C_fnthm, x, n, g, w, narm, drop, ret, nthreads)
# }
#
# fnthlC <- function(x, n = 0.5, g = NULL, w = NULL, narm = TRUE, drop = TRUE, ret = 1L, nthreads = 1L) {
#   .Call(C_fnthl, x, n, g, w, narm, drop, ret, nthreads)
# }


fquantile <- function(x, probs = c(0, 0.25, 0.5, 0.75, 1), w = NULL,
                      o = if(length(x) > 1e5L && length(probs) > log(length(x))) radixorder(x) else NULL,
                      na.rm = .op[["na.rm"]], type = 7L, names = TRUE,
                      check.o = is.null(attr(o, "sorted")))
  .Call(C_fquantile, x, probs, w, o, na.rm, type, names, check.o)

.quantile <- function(x, probs = c(0, 0.25, 0.5, 0.75, 1), w = NULL,
                      o = NULL, na.rm = TRUE, type = 7L, names = FALSE, check.o = FALSE)
  .Call(C_fquantile, x, probs, w, o, na.rm, type, names, check.o)

fscaleCpp <- function(x, ng = 0L, g = 0L, w = NULL, narm = TRUE, set_mean = 0, set_sd = 1) {
    .Call(Cpp_fscale, x, ng, g, w, narm, set_mean, set_sd)
}

fscalemCpp <- function(x, ng = 0L, g = 0L, w = NULL, narm = TRUE, set_mean = 0, set_sd = 1) {
    .Call(Cpp_fscalem, x, ng, g, w, narm, set_mean, set_sd)
}

fscalelCpp <- function(x, ng = 0L, g = 0L, w = NULL, narm = TRUE, set_mean = 0, set_sd = 1) {
    .Call(Cpp_fscalel, x, ng, g, w, narm, set_mean, set_sd)
}

fsumC <- function(x, ng = 0L, g = 0L, w = NULL, narm = TRUE, fill = FALSE, nthreads = 1L) {
    .Call(C_fsum, x, ng, g, w, narm, fill, nthreads)
}

fsummCcc <- function(x, w = NULL, drop = TRUE) {
  .Call(C_fsumm, x, 0L, 0L, w, FALSE, FALSE, drop, 1L)
}

fvarsdCpp <- function(x, ng = 0L, g = 0L, gs = NULL, w = NULL, narm = TRUE, stable_algo = TRUE, sd = TRUE) {
    .Call(Cpp_fvarsd, x, ng, g, gs, w, narm, stable_algo, sd)
}

fvarsdmCpp <- function(x, ng = 0L, g = 0L, gs = NULL, w = NULL, narm = TRUE, stable_algo = TRUE, sd = TRUE, drop = TRUE) {
    .Call(Cpp_fvarsdm, x, ng, g, gs, w, narm, stable_algo, sd, drop)
}

fvarsdlCpp <- function(x, ng = 0L, g = 0L, gs = NULL, w = NULL, narm = TRUE, stable_algo = TRUE, sd = TRUE, drop = TRUE) {
    .Call(Cpp_fvarsdl, x, ng, g, gs, w, narm, stable_algo, sd, drop)
}

mrtl <- function(X, names = FALSE, return = "list") {
  switch(return,
         list = .Call(Cpp_mrtl, X, names, 0L),
         data.frame = .Call(Cpp_mrtl, X, names, 1L),
         data.table = alc(.Call(Cpp_mrtl, X, names, 2L)),
         stop("Unknown return option!"))
}

mctl <- function(X, names = FALSE, return = "list") {
  switch(return,
         list = .Call(Cpp_mctl, X, names, 0L),
         data.frame = .Call(Cpp_mctl, X, names, 1L),
         data.table = alc(.Call(Cpp_mctl, X, names, 2L)),
         stop("Unknown return option!"))
}

psmatCpp <- function(x, g, t = NULL, transpose = FALSE, fill = NULL) {
    .Call(Cpp_psmat, x, g, t, transpose, fill)
}

qFCpp <- function(x, ordered = TRUE, na_exclude = TRUE, keep_attr = TRUE, ret = 1L) {
  .Call(Cpp_qF, x, ordered, na_exclude, keep_attr, ret)
}

sortuniqueCpp <- function(x) {
    .Call(Cpp_sortunique, x)
}

fdroplevelsCpp <- function(x, check_NA = TRUE) {
  .Call(Cpp_fdroplevels, x, check_NA)
}


setAttributes <- function(x, a) .Call(C_setAttributes, x, a)


copyMostAttributes <- function(to, from) .Call(C_copyMostAttributes, to, from)


setattributes <- function(x, a) .Call(C_setattributes, x, a) # invisible()


duplAttributes <- function(x, y) .Call(C_duplAttributes, x, y)


# No longer needed...
# setattr <- function(x, a, v) {
#   invisible(.Call(C_setattr, x, a, v))
# }

# duplattributes <- function(x, y) {
#     invisible(.Call(C_duplattributes, x, y))
# }

# cond_duplAttributes <- function(x, y) {
#     .Call(C_cond_duplAttributes, x, y)
# }

# cond_duplattributes <- function(x, y) {
#     invisible(.Call(C_cond_duplattributes, x, y))
# }

seqid <- function(x, o = NULL, del = 1L, start = 1L, na.skip = FALSE, skip.seq = FALSE, check.o = TRUE) {
  .Call(Cpp_seqid, x, o, del, start, na.skip, skip.seq, check.o)
}

groupid <- function(x, o = NULL, start = 1L, na.skip = FALSE, check.o = TRUE) {
  .Call(Cpp_groupid, x, o, start, na.skip, check.o)
}

