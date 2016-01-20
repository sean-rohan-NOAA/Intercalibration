#include <TMB.hpp>

// Compare two gear types

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* Data and parameters */
  DATA_VECTOR(SizeClass);     // Length in each length class
  DATA_ARRAY(N);                // haul x size
  DATA_VECTOR(SweptArea);       // Sample area for each haul. length=nhaul
  DATA_FACTOR(group);           // length(group)=nhaul. Defines pair ids
  DATA_FACTOR(Gear);            // length(Gear)=nhaul. Defines the two gear types
  DATA_SCALAR(huge);            // "Infinite" standard deviation on log scale
  DATA_SCALAR(tiny);            // "Zero" standard deviation on log scale

  PARAMETER_ARRAY(logspectrum); // group x size. One spectrum for each pair id
  PARAMETER_ARRAY(residual);    // haul x size. One AR residual for each haul

  // Selectivity curve
  PARAMETER(logsel0);     // Reletive selectivity in very small fish
  PARAMETER(logsel1);     // Reletive selectivity in very large fish
  PARAMETER(logselL50);     // Position of transition of the selectivity curve
  PARAMETER(logselSpread);     // Width of transition of the selectivity curve
  
  PARAMETER(logsd);             // sd of RW increments for the size spectrum at station
  PARAMETER(phi);               // AR1(phi) contribution of residuals
  PARAMETER(logsdres);          // sd of residuals
  PARAMETER(alpha);             // Exponent on the effect of SweptArea

  // Transpose for access by row:
  // array<Type> tnugget=nugget.transpose();
  array<Type> tresidual=residual.transpose();
  array<Type> tlogspectrum=logspectrum.transpose();
  array<Type> tN=N.transpose();

  vector<Type> loggear(N.dim[1]);

  int nhaul=N.dim[0];
  int nsize=N.dim[1];
  int ngear=NLEVELS(Gear);
  Type ans=0;
  Type sd=exp(logsd);

  // "Huge" or "Strictly infinte" variance on first size group?
  ans -= dnorm(vector<Type>(logspectrum.col(0)),Type(0),huge,true).sum();

  // Random walk over size spectrum at each station
  for(int i=1;i<nsize;i++){
    ans -= dnorm(vector<Type>(logspectrum.col(i)-logspectrum.col(i-1)),Type(0),sd,true).sum();
  }

  // AR(1) residuals
  using namespace density;
  SCALE_t< AR1_t<N01<Type> > >  nldens=SCALE(AR1(phi),exp(logsdres));
  for(int i=0;i<nhaul;i++){
    ans+=nldens(tresidual.col(i));
  }

  // Determine log gear effects
  Type sel0 = exp(logsel0);
  Type sel1 = exp(logsel1);
  Type selL50 = exp(logselL50);
  Type selSpread = exp(logselSpread);

  loggear = log(sel0 + (sel1-sel0)/(Type(1)+exp( (selL50 - SizeClass)/logselSpread)))*Type(0.5);

  ADREPORT(loggear);

  // Element-wise
  // for(int j=1;j<nsize;++j)
  //   {
  //     loggear(j) = 0.5*log(sel0 + (sel1-sel0)*1/(1+exp( (selL50 - SizeClass(j))/logselSpread)));
  //   }

  // Random walk prior on gear effect
  // ans-=dnorm(loggear(i,0),Type(0),huge,true);
  // for(int j=1;j<nsize;j++)
  // ans -= dnorm(loggear(j)-loggear(j-1),Type(0),sdGearRW,true);

  // Add data
  vector<Type> logintensity(nsize);
  for(int i=0;i<nhaul;i++)
    {
      logintensity=
	tlogspectrum.col(group[i])
	+tresidual.col(i) 
	+alpha*log(SweptArea(i));

      if(Gear(i)==1)
	{
	  logintensity += loggear;
	}
      else
	logintensity -= loggear;
      
      ans-=dpois(vector<Type>(tN.col(i)),exp(logintensity),true).sum();
    }

  return ans;
}

