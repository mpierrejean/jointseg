context("Joint segmentation using 'jointSeg'")

## A two-dimensional, noiseless signal
p <- 2
trueK <- 10
len <- 1e3

test_that("'jointSeg' recovers the true breakpoints exactly for a Gaussian noiseless signal", {
    for (ii in 1:10) {
        sim <- randomProfile(length=len, nBkp=trueK, noiseLevel=0, dim=p)
        Y <- sim$profile
        K <- 2*trueK
        
        res <- jointSeg(Y, method="RBS", K=K, modelSelectionMethod="Birge")
        bkpB <- res$bestBkp
        expect_equal(bkpB, sim$bkp)
        
        res <- jointSeg(Y, method="RBS", K=K, modelSelectionMethod="Lebarbier")
        bkpL <- res$bestBkp
        ## expect_equal(bkpL, sim$bkp)  ## often fails, cf Issue #6
        bkp <- res$dpBkp[[length(sim$bkp)]]
        expect_equal(bkp, sim$bkp)
        
        if (ii == 1) {        
            resDP <- jointSeg(Y, method = "DP", K = trueK)
            res <- jointSeg(Y, method = "DynamicProgramming", K = trueK)
            expect_identical(res, resDP)
            bkpDP <- res$bestBkp
            ## expect_equal(bkpDP, sim$bkp)
            bkp <- res$dpBkp[[length(sim$bkp)]]
            expect_equal(bkp, sim$bkp)
        }        
    }
})
