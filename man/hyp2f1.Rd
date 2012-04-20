\encoding{utf-8}
\name{hyp2f1}
\alias{hyp2f1}
\title{Compute the Gaussian hypergeometric function with complex arguments}
\usage{
  hyp2f1(a, b, c, z,
    algorithm = c("michel.stoitsov", "forrey"))
}
\arguments{
  \item{a}{complex parameter}

  \item{b}{complex parameter}

  \item{c}{complex parameter}

  \item{z}{complex variable}

  \item{algorithm}{either \dQuote{michel.stoitsov}
  (default) or \dQuote{forrey} (see the details)}
}
\value{
  The complex value of the Gaussian hypergeometric
  function.
}
\description{
  Two different algorithms can be used.
}
\details{
  The first, default, algorithm uses Fortran code in
  \dQuote{hyp_2F1.f90} from N. L. J.  Michel and M. V.
  Stoitsov, which is available at
  \url{http://cpc.cs.qub.ac.uk/summaries/AEAE}. The
  corresponding background reference is N. L. J. Michel and
  M. V. Stoitsov (2008): Fast computation of the Gauss
  hypergeometric function with all its parameters complex
  with application to the Pöschl-Teller-Ginocchio potential
  wave functions, Computer Physics Communications
  178:535-551.

  The second algorithm uses Fortran code in \dQuote{cyp.f}
  from R. C. Forrey is used which is available at
  \url{http://physics.bk.psu.edu/codes/chyp.f}. The
  corresponding background reference is R. C. Forrey
  (1997): Computing the hypergeometric function, Journal of
  Computational Physics 137:79-100.
}
\examples{
## library(appell)

## compare the results of both algorithms
## for random test data.

## todo: add better tests trying to replicate published results? 

nTest <- 100L
set.seed(38)

a <- complex(real=rnorm(nTest),
             imaginary=rnorm(nTest))
b <- complex(real=rnorm(nTest),
             imaginary=rnorm(nTest))
c <- complex(real=rnorm(nTest),
             imaginary=rnorm(nTest))
z <- complex(real=rnorm(nTest),
             imaginary=rnorm(nTest))

tableHyp2f1 <- matrix(nrow=nTest,
                      ncol=2L,
                      dimnames=
                      list(NULL,
                           c("forrey", "michel.stoitsov")))

for(i in seq_len(nTest))
{
    tableHyp2f1[i, "forrey"] <- hyp2f1(a[i], b[i], c[i], z[i],
                                       algorithm="forrey")
    tableHyp2f1[i, "michel.stoitsov"] <- hyp2f1(a[i], b[i], c[i], z[i],
                                                algorithm="michel.stoitsov")
}

tableHyp2f1

abs(tableHyp2f1[, "forrey"] - tableHyp2f1[, "michel.stoitsov"])
## so very small differences,
## at least in this range of function parameters.
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{math}

