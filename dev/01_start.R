# Building a Prod-Ready, Robust Shiny Application.
# 
# README: each step of the dev files is optional, and you don't have to 
# fill every dev scripts before getting started. 
# 01_start.R should be filled at start. 
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
# 
# 
########################################
#### CURRENT FILE: ON START SCRIPT #####
########################################

## Fill the DESCRIPTION ----
## Add meta data about your application
## 
## /!\ Note: if you want to change the name of your app during development, 
## either re-run this function, call golem::set_golem_name(), or don't forget
## to change the name in the app_sys() function in app_config.R /!\
## 
golem::fill_desc(
  pkg_name = "WimtrapWeb", # The Name of the package containing the App 
  pkg_title = "Plant 'Cis-Acting Regulatory Elements' Prediction Tools", # The Title of the package containing the App 
  pkg_description = "Application implementing the Wimtrap method to predict in different plant species the Cis-Acting Regulatory Elements
  and gene targets of transcription factors. The method inttegrates condition-independent and condition-depentdent genomic features at
  location of motif occurrences and makes use of TF-ChIP data to build predictive models based on supervised machine learning. Condition
  -independent genomic features are typicall the results of pattern-matching, the location on the genome and the conservation of the sequence.
  Condition-dependent genomic features include transcription factor footprints and variables related to the chromatin state.", # The Description of the package containing the App 
  author_first_name = "Quentin", # Your First Name
  author_last_name = "Rivière", # Your Last Name
  author_email = "qri@hotmail.be", # Your Email
  repo_url = "https://github.com/RiviereQuentin/WimtrapWeb" # The URL of the GitHub Repo (optional) 
)     

## Set {golem} options ----
golem::set_golem_options()

## Create Common Files ----
## See ?usethis for more information
usethis::use_gpl_license(version = 2)  # You can set another license here
usethis::use_readme_rmd( open = FALSE )
#usethis::use_code_of_conduct()
#usethis::use_lifecycle_badge( "Experimental" )
#usethis::use_news_md( open = FALSE )

## Use git ----
usethis::use_git()

## Init Testing Infrastructure ----
## Create a template for tests
golem::use_recommended_tests()

## Use Recommended Packages ----
golem::use_recommended_deps()

## Favicon ----
# If you want to change the favicon (default is golem's one)
#golem::use_favicon("C:/Users/qri/Documents/Carepat/Reviews/PlantCarepat_Logo3.png") # path = "path/to/ico". Can be an online file. 
golem::remove_favicon()

## Add helper functions ----
golem::use_utils_ui()
golem::use_utils_server()

# You're now set! ----

# go to dev/02_dev.R
rstudioapi::navigateToFile( "dev/02_dev.R" )

