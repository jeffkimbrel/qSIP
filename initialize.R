devtools::load_all()
devtools::check()
devtools::document()

devtools::install("../qSIP/")
packageVersion("qSIP")

# build package down site
#usethis::use_pkgdown()
usethis::use_pkgdown_github_pages()
pkgdown::build_site()


