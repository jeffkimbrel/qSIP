devtools::load_all()
devtools::check()
devtools::document()



# build package down site
#pkgdown::clean_site()
#usethis::use_pkgdown()
#usethis::use_pkgdown_github_pages()
#usethis::use_github_action("pkgdown")
pkgdown::build_site() ## change version in description first



devtools::install("../qSIP/")
packageVersion("qSIP")
