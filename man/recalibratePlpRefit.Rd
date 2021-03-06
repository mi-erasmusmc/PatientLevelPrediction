% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Recalibration.R
\name{recalibratePlpRefit}
\alias{recalibratePlpRefit}
\title{recalibratePlpRefit}
\usage{
recalibratePlpRefit(plpModel, newPopulation, newData, testFraction = 0.25)
}
\arguments{
\item{plpModel}{The trained plpModel (runPlp$model)}

\item{newPopulation}{The population created using createStudyPopulation() who will have their risks predicted}

\item{newData}{An object of type \code{plpData} - the patient level prediction
data extracted from the CDM.}

\item{testFraction}{Fraction of data to used for internal validation}
}
\value{
An object of class \code{runPlp} that is recalibrated on the new data
}
\description{
Train various models using a default parameter gird search or user specified parameters
}
\details{
The user can define the machine learning model to train (regularised logistic regression, random forest,
gradient boosting machine, neural network and )
}
