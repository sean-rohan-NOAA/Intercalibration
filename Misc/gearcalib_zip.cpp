#include <TMB.hpp>

// Random Walk density ('huge' is optional)
template<class Type>
  Type RW_logdens(vector<Type> x, Type sd, int order=1){
    for(int i=0; i<order; i++) x = diff(x);
    return dnorm(x, Type(0), sd, true).sum();
  }
template<class Type>
  Type RW_logdens(vector<Type> x, Type sd, Type huge, int order=1){
    Type ans = 0;
    for(int i=0; i<order; i++) ans += dnorm(x(i), Type(0), huge, true);
    ans += RW_logdens(x, sd, order);
    return ans;
  }

// Compare two gear types

template<class Type>
  Type objective_function<Type>::operator() ()
{
  /* Data and parameters */
  DATA_INTEGER(model_type); // 0=Poisson, 1=NegBin, 2=ZIP
  DATA_ARRAY(N);                // haul x size
  DATA_VECTOR(SweptArea);       // Sample area for each haul. length=nhaul
  DATA_FACTOR(group);           // length(group)=nhaul. Defines pair ids
  DATA_FACTOR(Gear);            // length(Gear)=nhaul. Defines the two gear types
  DATA_SCALAR(huge);            // "Infinite" standard deviation on log scale
  DATA_SCALAR(tiny);            // "Zero" standard deviation on log scale
  DATA_IVECTOR(rw_order);       // Order of logspectrum and loggear RW
  PARAMETER_ARRAY(logspectrum); // group x size. One spectrum for each pair id
  PARAMETER_ARRAY(nugget);      // haul x size. One nugget for each haul
  PARAMETER_ARRAY(residual);    // haul x size. One AR residual for each haul
  PARAMETER_VECTOR(loggear);     // 1 x size
  PARAMETER(logsd);             // sd of RW increments
  PARAMETER(phi);               // AR1(phi) contribution of residuals
  PARAMETER(logsdnug);          // sd of nugget effect
  PARAMETER(logsdres);          // sd of residuals
  PARAMETER(logsdGearRW);       // sd of increments in gear effect
  PARAMETER(logalpha);             // Exponent on the effect of SweptArea
  PARAMETER(logtheta);     // Negative binomial model log theta
  PARAMETER(logitpi);          // Zero-inflated Poisson logit pi
  
  // Transpose for access by row:
    array<Type> tnugget=nugget.transpose();
  array<Type> tresidual=residual.transpose();
  array<Type> tlogspectrum=logspectrum.transpose();
  array<Type> tN=N.transpose();
  
  int nhaul=N.dim[0];
  int nsize=N.dim[1];
  int ngear=NLEVELS(Gear);
  Type ans=0;
  Type sd=exp(logsd);
  Type sdGearRW=exp(logsdGearRW);
  Type sdnug=exp(logsdnug);
  Type alpha=exp(logalpha);
  Type theta=exp(logtheta);
  Type pi = 1 / (1 + exp(-logit_pi)); // Inverse logit to get pi in (0, 1)
  
  // Random walk over size spectrum at each station
  for(int i=0; i<tlogspectrum.cols(); i++){
    ans -= RW_logdens(vector<Type>(tlogspectrum.col(i)), sd, huge, rw_order(0));
  }
  
  // AR(1) residuals
  using namespace density;
  SCALE_t< AR1_t<N01<Type> > >  nldens=SCALE(AR1(phi),exp(logsdres));
  for(int i=0; i<nhaul; i++){
    ans += nldens(tresidual.col(i));
  }
  
  // White noise nugget effect
  for(int j=0;j<nsize;j++)
    ans -= dnorm(vector<Type>(nugget.col(j)),Type(0),sdnug,true).sum();
  
  // Random walk prior on gear effect
  ans -= RW_logdens(loggear, sdGearRW, huge, rw_order(1));
  
  // Add data
  vector<Type> logintensity(nsize);
  for(int i=0;i<nhaul;i++)
  {
    logintensity=
      tlogspectrum.col(group[i])
    +tresidual.col(i) 
    +alpha*log(SweptArea(i))
    +tnugget.col(i);
    if(Gear(i)==1)
    {
      logintensity += loggear;
    }
    else
      logintensity -= loggear;
      
      vector<Type> mu = exp(logintensity);
      
      for (int j = 0; j < nsize; j++) {
        Type y = tN(i, j);
        
        if (model_type ==  1) {
          // Poisson
          ans -= dpois(y, mu(j), true);
          
        } else if (model_type == 2) {
          // Negative Binomial
          Type var = mu(j) + mu(j) * mu(j) / theta;
          ans -= dnbinom2(y, mu(j), var, true);
          
        } else if (model_type == 3) {
          // Zero-inflated Poisson
          if (y == 0) {
            ans -= log(pi + (1 - pi) * exp(dpois(Type(0), mu(j), true)));
          } else {
            ans -= log(1 - pi) + dpois(y, mu(j), true);
          }
          
        } else {
          // Safety fallback (should never be hit if model_type is validated externally)
          error("Invalid model_type specified.");
        }
      }
  }
  
  return ans;
  }