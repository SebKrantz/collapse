# Make package website:
# usethis::use_pkgdown() # overwrites _pkgdown.yml
library(pkgdown) # Built with pkgdown 1.6.1
# update .pkgdown.yml
init_site()
build_home(preview = FALSE)
# in index.html, remove misc/figures/.
topics <- sort(setdiff(unlist(lapply(tools::Rd_db("collapse"),
                                    tools:::.Rd_get_metadata, "name"), use.names = FALSE),
                      c("collapse-documentation","A0-collapse-documentation","collapse-depreciated","collapse-renamed"))) # "collapse-package"
# TODO: Update plm method documentation !!!
options(max.print = 100L)
build_reference(examples = TRUE, topics = topics) # "collapse-package"
Sys.setenv(NCRAN = "TRUE")
Sys.setenv(RUNBENCH = "TRUE")
# build_articles(lazy = TRUE) # lazy = FALSE # Still do with NCRAN = TRUE
build_articles_index()
build_article("collapse_intro")
build_article("collapse_and_plm")
build_article("collapse_and_data.table")
build_article("collapse_and_sf")
# Replace all A0-collapse-documentation.html with index.html !!
# Also replace A1-fast-statistical-functions.html with fast-statistical-functions.html etc.
# also replace / remove %20 (empty spaces) (nah...)
# Replce all <h1>Reference</h1> with <h1>Documentation & Overview</h1>
# Replce all <h1>Articles</h1> with <h1>Vignettes / Articles</h1>
# Replace &amp;lt;- with &lt;- and %&amp;gt;% with %&gt;% and
# <span class='kw'>&amp;</span><span class='no'>lt</span>;- with <span class='kw'>&lt;-</span>
build_news()
# Add favicon:
# replace <head> by <head><link rel="shortcut icon" href="https://sebkrantz.github.io/collapse/favicon.ico" type="image/x-icon">
# or <head><meta

#     <a href = "collapse1.7.digest.html"><strong><em>collapse</em> 1.7 NEWS Digest</strong></a>

# replace "https://twitter.com/collapse_R" target="_blank",
# href="https://github.com/SebKrantz/collapse" target="_blank" and
# href="https://sebkrantz.github.io/Rblog/" target="_blank"  (target = "_blank")

# in index.html, replace with
# <li id="bloglink">
#  <a href="https://sebkrantz.github.io/Rblog/" target="_blank">Blog</a>
# </li>

# in pkgdown.js, add line
# $("#bloglink").removeClass("active");   below     menu_anchor.closest("li.dropdown").addClass("active");

# Finally: Fix favicons by replacing "fas fa- with "fa fa-
preview_site()



# ?build_home
