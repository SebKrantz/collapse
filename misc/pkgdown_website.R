# Make package website:
# usethis::use_pkgdown() # overwrites _pkgdown.yml
library(pkgdown)
# init_site()
build_home(preview = FALSE)
# in index.html, remove misc/figures.
topics <- sort(setdiff(unlist(lapply(tools::Rd_db("collapse"),
                                    tools:::.Rd_get_metadata, "name"), use.names = FALSE),
                      c("collapse-documentation","A0-collapse-documentation","collapse-depreciated"))) # "collapse-package"
build_reference(examples = TRUE, topics = topics) # "collapse-package"
# build_articles(lazy = FALSE) # Still do with NCRAN = TRUE
# Replace all collapse-documentation.html with index.html !!
# Replce all <h1>Reference</h1> with <h1>Documentation & Overview</h1>
# Replce all <h1>Articles</h1> with <h1>Vignettes / Articles</h1>
# Replace &amp;lt;- with &lt;- and %&amp;gt;% with %&gt;% and
# <span class='kw'>&amp;</span><span class='no'>lt</span>;- with <span class='kw'>&lt;-</span>
build_news()
# Add favicon:
# replace <head> by <head><link rel="shortcut icon" href="https://sebkrantz.github.io/collapse/favicon.ico" type="image/x-icon">
preview_site()

# ?build_home
