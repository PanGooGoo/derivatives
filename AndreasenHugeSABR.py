from _library import IRAnalytics
import numpy
from scipy.stats import norm

class SABRAndreasenHuge():
    
    def __init__(self, alpha, beta, rho, nu, expiry, strikes, fwd, strike_step = 0):
        self.alpha = alpha
        self.beta = beta
        self.rho = rho
        self.nu = nu
        self.expiry = expiry
        self.strikes = strikes  
        self.numb_strikes = len(self.strikes)
        self.fwd = fwd
        self.strike_step = strike_step
        self.y = 0
        self.J = 0
        self.x = 0
        self.expansionNormVol = 0
        self.expansionLogNormVol = 0
        self.expansionLocalVol = 0
        self.lv_call_px = numpy.empty(self.numb_strikes)
        self.norm_call_px = numpy.empty(self.numb_strikes)
        self.analytics = IRAnalytics.IRAnalytics()
        self.computeExpansionTerms()
                    
    def computeExpansionTerms(self):
        try:
            self.y = (numpy.power(self.fwd,(1.0 - self.beta)) - numpy.power(self.strikes,(1.0 - self.beta)))/(self.alpha * (1.0 - self.beta))
            self.J = numpy.sqrt(1.0 + self.nu * self.nu * self.y * self.y - 2.0 * self.rho * self.nu * self.y)
            self.x = (1.0 / self.nu) * numpy.log((self.J + self.nu * self.y - self.rho) / (1.0 - self.rho))
        except Exception as e:
            print 'Error in computeExpansionTerms: %s' %(e)
            
    def computeExpansionNormalVol(self):
        try:
            ret = numpy.empty(self.numb_strikes)
            i = 0
            for strike in self.strikes:
                if strike == self.fwd:
                    ret[i] = numpy.power(self.fwd, self.beta) * self.alpha
                else:
                    ret[i] = (self.fwd - strike) / self.x[i]
                i += 1
            self.expansionNormVol = ret
        except Exception as e:
            print 'Error in computeExpansionNormalVol: %s' %(e)
            
    def getExpansionNormalVol(self):
        return self.expansionNormVol
    
    def computeExpansionNormalVolPrice(self):
        try:
            self.computeExpansionNormalVol()
            i = 0
            for vol in self.expansionNormVol:
                self.norm_call_px[i] = self.analytics.NormalBlackScholes(self.fwd, self.strikes[i], vol, "Call", self.expiry)
                i += 1
        except Exception as e:
            print 'Error in computeExpansionNormalVolPrice: %s' %(e)

    def getExpansionNormalVolPrice(self):
        return self.norm_call_px
                        
    def computeExpansionLogNormalVol(self):
        try:
            ret = numpy.empty(self.numb_strikes)
            i = 0
            for strike in self.strikes:
                if strike == self.fwd:
                    ret[i] = numpy.power(self.fwd, (self.beta - 1.0)) * self.alpha
                else:
                    ret[i] = numpy.log(self.fwd / strike) / self.x[i]
                i += 1
            self.expansionLogNormVol = ret
        except Exception as e:
            print 'Error in computeExpansionLogNormalVol: %s' %(e)
            
    def getExpansionLogNormalVol(self):
        return self.expansionLogNormVol
                            
    def computeExpansionLocalVol(self):
        try:
            self.expansionLocalVol = self.J * numpy.power(self.strikes, self.beta) * self.alpha 
        except Exception as e:
            print 'Error in computeExpansionLocalVol: %s' %(e)

    def getExpansionLocalVol(self):
        return self.expansionLocalVol
    
    def computeExpansionLocalVolPrice(self, use_lv_adjustment=False):
        try:
            self.computeExpansionLocalVol()
            lv_adj = 1.0
            if use_lv_adjustment:
                lv_adj = self.computeLocalVolAdjustment()
            theta_sqrd = numpy.power(self.getExpansionLocalVol(), 2.0) 
            theta_sqrd = numpy.multiply(theta_sqrd, lv_adj)
            strike_step_sqrd = self.strike_step*self.strike_step

            # Implementation based on Volatility Interpolation - Andreasen, Huge
            A = numpy.zeros((self.numb_strikes,self.numb_strikes))
            i = 1
            j = 0
            while i < self.numb_strikes - 1:
                A[i,j] = -0.5*self.expiry/strike_step_sqrd*theta_sqrd[i]
                A[i,j+1] = 1.0+self.expiry/strike_step_sqrd*theta_sqrd[i]
                A[i,j+2] = -0.5*self.expiry/strike_step_sqrd*theta_sqrd[i]
                i += 1
                j += 1
            A[0,0] = 1
            A[self.numb_strikes-1,self.numb_strikes-1] = 1
            self.computeExpansionNormalVol()
            normal_vols = self.expansionNormVol
            B = numpy.maximum(self.fwd - self.strikes,0)
            B[0] = self.analytics.NormalBlackScholes(self.fwd, self.strikes[0], normal_vols[0], "Call", self.expiry)
            B[self.numb_strikes-1] = self.analytics.NormalBlackScholes(self.fwd, self.strikes[self.numb_strikes-1], normal_vols[self.numb_strikes-1], "Call", self.expiry)
            self.lv_call_px = numpy.linalg.solve(A, B)

        except Exception as e:
            print 'Error in computeExpansionLocalVolPrice: %s' %(e)
            
    def getExpansionLocalVolPrice(self):
        return self.lv_call_px
    
    def computeLocalVolAdjustment(self):
        try:
            xi = abs(self.x)/numpy.sqrt(self.expiry)
            lv_adjustment = 2.0 * (1.0 - xi*norm.cdf(-xi)/norm.pdf(xi))
            return lv_adjustment

        except Exception as e:
            print 'Error in computeLocalVolAdjustment: %s' %(e)
            
    def computePDF(self, prices):
        try:
            p = numpy.zeros(self.numb_strikes)
            i = 1
            while ((i > 0) and (i < self.numb_strikes-1)):
                p[i] = (prices[i-1]-2.0*prices[i]+prices[i+1])/(self.strike_step*self.strike_step)
                i += 1
            return p
            
        except Exception as e:
            print 'Error in computePDF: %s' %(e)
                        
if __name__ == '__main__':
    '''
#     MC vs. AH test
#     the results for this test are reported in the ZAB.xls file
    alpha = 0.02
    beta = 0.5
    rho = 0.4
    nu = 0.5
    expiry = 19.5
    fwd = 0.025
    strike_step = 0.0025
    strikes = numpy.arange(0.0, 0.0525, strike_step)
    '''
#     AH Paper Parameters
#     the results for this test are reported in the ZAB.xls file
    alpha = 0.0873
    beta = 0.7
    rho = -0.48
    nu = 0.47
    expiry =10.0
    fwd = 0.0325
    strike_step = 0.0025
    strikes = numpy.arange(0.0, 0.1025, strike_step)

    print 'Strikes:'
    print strikes
    sabrah = SABRAndreasenHuge(alpha, beta, rho, nu, expiry, strikes, fwd, strike_step)
    sabrah.computeExpansionLocalVolPrice(use_lv_adjustment=True)
    print 'Px from Local Vol (matrix form):'
    px_from_lv = sabrah.getExpansionLocalVolPrice()
    print px_from_lv
    print 'Px from Local Vol:'
    for item in px_from_lv:
        print item
    print 'PDF from Local Vol (matrix form):'
    lv_pdf = sabrah.computePDF(px_from_lv)
    print lv_pdf
    print 'PDF from Local Vol:'
    for item in lv_pdf:
        print item
#     print 'Expansion Px from Normal Vol:'
#     sabrah.computeExpansionNormalVolPrice()
#     for item in sabrah.getExpansionNormalVolPrice():
#         print item