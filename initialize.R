devtools::load_all()
devtools::check()
devtools::document()

devtools::install("../qSIP/")
packageVersion("qSIP")

# build package down site
usethis::use_pkgdown()
pkgdown::build_site()
