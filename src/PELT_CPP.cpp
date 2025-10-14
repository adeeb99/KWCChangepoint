// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::uvec mySeq(unsigned int s, unsigned int e){

  unsigned int len=e-s+1;
  arma::uvec pos(len);

  for(unsigned int i=0;i<len;i++){
    pos[i]=s+i;
  }

  return pos;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double cost_cpp(unsigned int s, unsigned int e,arma::rowvec data){

  double len=e-s+1;
  double sum=0;
  double val;

  for(unsigned int i=s;i<e+1;i++){
    sum=sum+ data[i];

  }

  val=-(pow(sum,2))/len;

  return val;
}

// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// double cost_sch(unsigned int s, unsigned int e,arma::rowvec data,double N){
//
//   double len=e-s+1;
//   double sum=0;
//   double mean=0;
//   double val;
//
//   for(unsigned int i=s;i<e+1;i++){
//     sum=sum+ data[i];
//
//   }
//
//   mean=sum/len;
//
//   sum=0;
//
//   for(unsigned int i=s;i<e+1;i++){
//     sum=sum+pow(data[i]-mean,2);
//
//   }
//
//   sum=sum/len;
//
//   val=-sum;
//
//   return val;
// }

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::uvec updateR(arma::uvec R1, arma::rowvec obj, double F_tau_star, double beta,double K,double tau_star){



  unsigned int counter=1;

  for(unsigned int i=0;i<R1.size();i++){


    if((obj[i]-beta+K)<=F_tau_star){

      counter++;
    }
  }

  arma::uvec newR(counter);

  newR[counter-1]=tau_star;
  counter=0;

  for(unsigned int i=0;i<R1.size();i++){
    if((obj[i]-beta+K)<=F_tau_star){
      newR[counter]=R1[i];
      counter++;
    }
  }

  return newR;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat PELT(arma::rowvec data, int N,double beta,double K=0) {

  // of a given size (filled with 0)
  arma::rowvec F_vals=zeros(1,N+1);

  double va=(pow(N,2)-1)/12;
//  double mu_N=va*3*(N+1)/N;



  F_vals[0]=-beta;

 // arma::uvec pos;

  arma::mat cp=zeros(N+1,N+1);
  // Rcout<<"CP "<< cp<< std::endl;
  double C;
  unsigned int tau_1;
  uword ob_min;

  arma::uvec R1(1);
  arma::rowvec obj;

  R1[0]=0;
 // Rcout << "The value F is " << F_vals << std::endl;
  for(int tau_star = 1; tau_star <N+1; tau_star++) {

    obj=zeros(1,R1.size());
    for(unsigned int i=0; i<R1.size(); i++){

          // if(!sch){
            C=cost_cpp(R1[i],tau_star-1,data);
            obj[i]=F_vals[R1[i]]+beta+C/va;
          // }
          // else
          // {
          //   C=cost_sch(R1[i],tau_star-1,data,N);
          //   obj[i]=F_vals[R1[i]]+beta+C;
          // }

    }



    ob_min=obj.index_min();

    tau_1=R1[ob_min];
    F_vals[tau_star]=obj[ob_min];



    //starts at 0,
    //print()
    cp.row(tau_star)=cp.row(tau_1);
    cp(tau_star,tau_1)=1;




    R1=updateR(R1,obj,F_vals[tau_star],beta,K,tau_star*1.0);
//    Rcout << "The value tau_1 is " << tau_1 << std::endl;
 //    Rcout << "cp " <<cp << std::endl;
  //   Rcout << "The value F_val is " << F_vals << std::endl;
  // // Rcout << "The value obj is " << obj << std::endl;
  // Rcout << "The value R1 is " << R1.t() << std::endl;
  // Rcout << "The value tau_star is " << tau_star << std::endl;
   // Rcout << "The value tau_1 is " << tau_1 << std::endl;
  }
 // Rcout << "The value F_val is " << F_vals[N] << std::endl;
  return cp.row(N);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int getOnes(arma::mat vec){
  int m=0;
  for(unsigned int i=0;i<vec.n_cols;i++){
    if(vec(0,i)==1){
      m++;
    }
  }
  return m;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat PELT_NL(arma::rowvec data, int N,double c2=1,int m0=1){

  arma::mat cp;
  double curr_beta=c2/m0;
  int m;

  cp=PELT(data,N,curr_beta,0);
  m=getOnes(cp);

  while(m0!=m){
    m0=m;
    curr_beta=c2/m;

    cp=PELT(data,N,curr_beta,0);
    m=getOnes(cp);


  }
  return cp;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R

set.seed(440)
test<-c(rnorm(50),rnorm(50,5),rnorm(50,3))
ranks<-rank(test)

ranks=1:20
which(PELT(ranks,length(ranks),1,0)==1)-1
##cost_cpp(0,10,ranks)
PELT_T<-function(ranks,beta,K=0){
  which(PELT(ranks,length(ranks),beta,0)==1)-1
}
set.seed(440)
test<-c(rnorm(50),rnorm(50,5),rnorm(50,3))
ranks<-rank(test)
which(PELT_NL(ranks,length(ranks),c2=10)==1)-1


PELT_T_NL<-function(ranks,c2,m0=1){
  which(PELT_NL(ranks,length(ranks),c2=c2,m0=m0)==1)-1
}

*/
