
BWCpp <- function(x, ng = 0L, g = 0L, gs = NULL, w = NULL, narm = TRUE, theta = 1, set_mean = 0, B = FALSE, fill = FALSE) {
    .Call(Cpp_BW, x, ng, g, gs, w, narm, theta, set_mean, B, fill)
}

BWmCpp <- function(x, ng = 0L, g = 0L, gs = NULL, w = NULL, narm = TRUE, theta = 1, set_mean = 0, B = FALSE, fill = FALSE) {
    .Call(Cpp_BWm, x, ng, g, gs, w, narm, theta, set_mean, B, fill)
}

BWlCpp <- function(x, ng = 0L, g = 0L, gs = NULL, w = NULL, narm = TRUE, theta = 1, set_mean = 0, B = FALSE, fill = FALSE) {
    .Call(Cpp_BWl, x, ng, g, gs, w, narm, theta, set_mean, B, fill)
}

TRACpp <- function(x, xAG, g = 0L, ret = 1L) {
    .Call(Cpp_TRA, x, xAG, g, ret)
}

TRAmCpp <- function(x, xAG, g = 0L, ret = 1L) {
    .Call(Cpp_TRAm, x, xAG, g, ret)
}

TRAlCpp <- function(x, xAG, g = 0L, ret = 1L) {
    .Call(Cpp_TRAl, x, xAG, g, ret)
}

fndistinctCpp <- function(x, ng = 0L, g = 0L, gs = NULL, narm = TRUE) {
    .Call(Cpp_fndistinct, x, ng, g, gs, narm)
}

fndistinctlCpp <- function(x, ng = 0L, g = 0L, gs = NULL, narm = TRUE, drop = TRUE) {
    .Call(Cpp_fndistinctl, x, ng, g, gs, narm, drop)
}

fndistinctmCpp <- function(x, ng = 0L, g = 0L, gs = NULL, narm = TRUE, drop = TRUE) {
    .Call(Cpp_fndistinctm, x, ng, g, gs, narm, drop)
}

pwnobsmCpp <- function(x) {
    .Call(Cpp_pwnobsm, x)
}

fnobsC <- function(x, ng = 0L, g = 0L) {
    .Call(C_fnobs, x, ng, g)
}

# fnobsmC <- function(x, ng = 0L, g = 0L, drop = TRUE) {
#     .Call(C_fnobsm, x, ng, g, drop)
# }

# fnobslC <- function(x, ng = 0L, g = 0L, drop = TRUE) {
#     .Call(C_fnobsl, x, ng, g, drop)
# }

varyingCpp <- function(x, ng = 0L, g = 0L, any_group = TRUE) {
    .Call(Cpp_varying, x, ng, g, any_group)
}

varyingmCpp <- function(x, ng = 0L, g = 0L, any_group = TRUE, drop = TRUE) {
    .Call(Cpp_varyingm, x, ng, g, any_group, drop)
}

varyinglCpp <- function(x, ng = 0L, g = 0L, any_group = TRUE, drop = TRUE) {
    .Call(Cpp_varyingl, x, ng, g, any_group, drop)
}

fbstatsCpp <- function(x, ext = FALSE, ng = 0L, g = 0L, npg = 0L, pg = 0L, w = NULL, array = TRUE, setn = TRUE, gn = NULL) {
    .Call(Cpp_fbstats, x, ext, ng, g, npg, pg, w, array, setn, gn)
}

fbstatsmCpp <- function(x, ext = FALSE, ng = 0L, g = 0L, npg = 0L, pg = 0L, w = NULL, array = TRUE, gn = NULL) {
    .Call(Cpp_fbstatsm, x, ext, ng, g, npg, pg, w, array, gn)
}

fbstatslCpp <- function(x, ext = FALSE, ng = 0L, g = 0L, npg = 0L, pg = 0L, w = NULL, array = TRUE, gn = NULL) {
    .Call(Cpp_fbstatsl, x, ext, ng, g, npg, pg, w, array, gn)
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

fmeanCpp <- function(x, ng = 0L, g = 0L, gs = NULL, w = NULL, narm = TRUE) {
    .Call(Cpp_fmean, x, ng, g, gs, w, narm)
}

fmeanmCpp <- function(x, ng = 0L, g = 0L, gs = NULL, w = NULL, narm = TRUE, drop = TRUE) {
    .Call(Cpp_fmeanm, x, ng, g, gs, w, narm, drop)
}

fmeanlCpp <- function(x, ng = 0L, g = 0L, gs = NULL, w = NULL, narm = TRUE, drop = TRUE) {
    .Call(Cpp_fmeanl, x, ng, g, gs, w, narm, drop)
}

fnthCpp <- function(x, n = 0.5, ng = 0L, g = 0L, gs = NULL, w = NULL, narm = TRUE, ret = 1L) {
  .Call(Cpp_fnth, x, n, ng, g, gs, w, narm, ret)
}

fnthmCpp <- function(x, n = 0.5, ng = 0L, g = 0L, gs = NULL, w = NULL, narm = TRUE, drop = TRUE, ret = 1L) {
  .Call(Cpp_fnthm, x, n, ng, g, gs, w, narm, drop, ret)
}

fnthlCpp <- function(x, n = 0.5, ng = 0L, g = 0L, gs = NULL, w = NULL, narm = TRUE, drop = TRUE, ret = 1L) {
  .Call(Cpp_fnthl, x, n, ng, g, gs, w, narm, drop, ret)
}

fmodeCpp <- function(x, ng = 0L, g = 0L, gs = NULL, w = NULL, narm = TRUE, ret = 0L) {
    .Call(Cpp_fmode, x, ng, g, gs, w, narm, ret)
}

fmodelCpp <- function(x, ng = 0L, g = 0L, gs = NULL, w = NULL, narm = TRUE, ret = 0L) {
    .Call(Cpp_fmodel, x, ng, g, gs, w, narm, ret)
}

fmodemCpp <- function(x, ng = 0L, g = 0L, gs = NULL, w = NULL, narm = TRUE, drop = TRUE, ret = 0L) {
    .Call(Cpp_fmodem, x, ng, g, gs, w, narm, drop, ret)
}

fprodCpp <- function(x, ng = 0L, g = 0L, w = NULL, narm = TRUE) {
    .Call(Cpp_fprod, x, ng, g, w, narm)
}

fprodmCpp <- function(x, ng = 0L, g = 0L, w = NULL, narm = TRUE, drop = TRUE) {
    .Call(Cpp_fprodm, x, ng, g, w, narm, drop)
}

fprodlCpp <- function(x, ng = 0L, g = 0L, w = NULL, narm = TRUE, drop = TRUE) {
    .Call(Cpp_fprodl, x, ng, g, w, narm, drop)
}

fscaleCpp <- function(x, ng = 0L, g = 0L, w = NULL, narm = TRUE, set_mean = 0, set_sd = 1) {
    .Call(Cpp_fscale, x, ng, g, w, narm, set_mean, set_sd)
}

fscalemCpp <- function(x, ng = 0L, g = 0L, w = NULL, narm = TRUE, set_mean = 0, set_sd = 1) {
    .Call(Cpp_fscalem, x, ng, g, w, narm, set_mean, set_sd)
}

fscalelCpp <- function(x, ng = 0L, g = 0L, w = NULL, narm = TRUE, set_mean = 0, set_sd = 1) {
    .Call(Cpp_fscalel, x, ng, g, w, narm, set_mean, set_sd)
}

fsumC <- function(x, ng = 0L, g = 0L, w = NULL, narm = TRUE) {
    .Call(C_fsum, x, ng, g, w, narm)
}

# fsummC <- function(x, ng = 0L, g = 0L, w = NULL, narm = TRUE, drop = TRUE) {
#     .Call(C_fsumm, x, ng, g, w, narm, drop)
# }

# fsumlC <- function(x, ng = 0L, g = 0L, w = NULL, narm = TRUE, drop = TRUE) {
#     .Call(C_fsuml, x, ng, g, w, narm, drop)
# }

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

psmatCpp <- function(x, g, t = NULL, transpose = FALSE) {
    .Call(Cpp_psmat, x, g, t, transpose)
}

qFCpp <- function(x, sort = TRUE, ordered = FALSE, na.exclude = TRUE, keep.attr = TRUE, ret = 1L) {
    .Call(Cpp_qF, x, sort, ordered, na.exclude, keep.attr, ret)
}

funiqueCpp <- function(x, sort = TRUE) {
    .Call(Cpp_funique, x, sort)
}

fdroplevelsCpp <- function(x, check_NA = TRUE) {
  .Call(Cpp_fdroplevels, x, check_NA)
}


setAttributes <- function(x, a) {
    .Call(C_setAttributes, x, a)
}

copyMostAttributes <- function(to, from) {
  .Call(C_copyMostAttributes, to, from)
}

setattributes <- function(x, a) {
    invisible(.Call(C_setattributes, x, a))
}

setattr <- function(x, a, v) {
    invisible(.Call(C_setattr, x, a, v))
}

duplAttributes <- function(x, y) {
    .Call(C_duplAttributes, x, y)
}

duplattributes <- function(x, y) {
    invisible(.Call(C_duplattributes, x, y))
}

cond_duplAttributes <- function(x, y) {
    .Call(C_cond_duplAttributes, x, y)
}

# cond_duplattributes <- function(x, y) {
#     invisible(.Call(C_cond_duplattributes, x, y))
# }

seqid <- function(x, o = NULL, del = 1L, start = 1L, na.skip = FALSE, skip.seq = FALSE, check.o = TRUE) {
  .Call(Cpp_seqid, x, o, del, start, na.skip, skip.seq, check.o)
}

groupid <- function(x, o = NULL, start = 1L, na.skip = FALSE, check.o = TRUE) {
  .Call(Cpp_groupid, x, o, start, na.skip, check.o)
}

