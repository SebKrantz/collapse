# Make package website:
# usethis::use_pkgdown() # overwrites _pkgdown.yml
library(pkgdown)
init_site()
build_home(preview = FALSE)
topics <- sort(setdiff(unlist(lapply(tools::Rd_db("collapse"),
                                    tools:::.Rd_get_metadata, "name"), use.names = FALSE),
                      c("collapse-documentation","A0-collapse-documentation","collapse-depreciated"))) # "collapse-package"
build_reference(examples = TRUE, topics = topics)
# build_articles(lazy = FALSE) # Still do with NCRAN = TRUE
# Replace all collapse-documentation.html with index.html !!
# Replce all <h1>Reference</h1> with <h1>Documentation & Overview</h1>
# Replce all <h1>Articles</h1> with <h1>Vignettes / Articles</h1>
build_news()
# in index.html, remove misc/figures.
preview_site()

# ?build_home
