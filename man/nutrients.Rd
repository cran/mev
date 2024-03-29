% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{nutrients}
\alias{nutrients}
\title{Nutrient data}
\format{
A data frame with 9544 rows and 38 variables:
\describe{
  \item{\code{prot}}{proteins (in grams)}
  \item{\code{carb}}{carbonhydrate (in gram)}
  \item{\code{sugr}}{total sugars (in gram)}
  \item{\code{fibe}}{dietary fibers (in grams)}
  \item{\code{tfat}}{total fat (in grams)}
  \item{\code{sfat}}{saturated fat (in grams)}
  \item{\code{mfat}}{monounsaturated fat (in grams)}
  \item{\code{pfat}}{polyunsaturated fat (in grams)}
  \item{\code{chol}}{cholesterol (in milligrams)}
  \item{\code{atoc}}{vitamin E as alpha-tocopherol (in milligrams)}
  \item{\code{ret}}{retinol (in micrograms)}
  \item{\code{vara}}{Vitamin A as retinol activity equivalents (in micrograms).}
  \item{\code{acar}}{alpha-carotene (in micrograms)}
  \item{\code{bcar}}{beta-carotene (in micrograms)}
  \item{\code{cryp}}{beta-cryptoxanthin (in micrograms)}
  \item{\code{lyco}}{lycopene (in micrograms)}
  \item{\code{lz}}{lutein and zeaxanthin (in micrograms).}
  \item{\code{vb1}}{thiamin (vitamin B1, in milligrams)}
  \item{\code{vb2}}{riboflavin (vitamin B2, in milligrams)}
  \item{\code{niac}}{niacin (in milligrams)}
  \item{\code{vb6}}{vitamin B5 (in milligrams)}
  \item{\code{fola}}{total folate (in micrograms)}
  \item{\code{fa}}{folic acid (in micrograms)}
  \item{\code{ff}}{food folate (in micrograms)}
  \item{\code{chl}}{total choline (in milligrams)}
  \item{\code{vb12}}{vitamin B12 (in micrograms)}
  \item{\code{vc}}{vitamin C (in milligrams)}
  \item{\code{vd}}{vitamin D (comprising D2 and D3, in micrograms)}
  \item{\code{vk}}{vitamin K (in micrograms)}
  \item{\code{calc}}{calcium (in milligrams)}
  \item{\code{phos}}{phosphorus (in milligrams)}
  \item{\code{magn}}{magnesium (in milligrams)}
  \item{\code{iron}}{iron (in milligrams)}
  \item{\code{zinc}}{zinc (in milligrams)}
  \item{\code{copp}}{copper (in milligrams)}
  \item{\code{sodi}}{sodium (in milligrams)}
  \item{\code{pota}}{potassium (in milligrams)}
  \item{\code{sele}}{selenium (in micrograms)}
}
}
\source{
National Center for Health Statistics, \url{https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/DR1TOT_I.XPT}
}
\usage{
nutrients
}
\description{
Interview component of survey 'What we eat in
America'. These are extracted from the 2015–2016 National Health and Nutrition Examination Survey (NHANES, \url{https://wwwn.cdc.gov/nchs/nhanes/Default.aspx}) report and consist of the total nutrients for all food and beverage intake ingested over a 24 hours period.
}
\details{
Note that the sample design oversampled specific population targets and that only respondants are provided. The website contains more information about sampling weights. There are multiple missing records.
}
\note{
These data are subject to a data user agreement, available at \url{https://www.cdc.gov/nchs/data_access/restrictions.htm}
}
\keyword{datasets}
