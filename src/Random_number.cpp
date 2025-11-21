// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <random>
#include <iostream>
#include <math.h>
//#include <ctime>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector random_unif_interval(int size, int max_val){

 // vector to store the generated random number
 NumericVector index(size);

  // initialize the random generator
  std::random_device seed;
  std::mt19937 rng(seed());

  // uncomment this line if std::random_device does not work properly
  //std::default_random_engine rng(time(0));

  std::uniform_int_distribution<> unif_c(1, max_val);
  for(int i=0; i<size; ++i){
    index[i] = unif_c(rng);
  }
  return(index);

}

// [[Rcpp::export]]
NumericVector random_unif_interval_seed(int size, int max_val, int input_seed){

  // create 20% more indices such that duplicated indices are less probable
  int new_size = size*1.2;
  NumericVector index_vector(new_size);

  // initialize the random generator
  std::mt19937 rng(input_seed);

  // define the std integer uniform distribution
  std::uniform_int_distribution<> unif_c(1, max_val);

  // generate "new_size" random number using the above distribution
  for(int i = 0 ; i < new_size ; i++){
    index_vector[i] = unif_c(rng);
  }

  // keep only unique numbers
  index_vector = unique(index_vector);
  // calculate how many unique number were generated
  int unique_size = index_vector.length();

  // if there amount of unique numbers is greater than the required size
  if(unique_size > size){
    return index_vector;
  }
  else{ // if there amount of unique numbers is NOT enough

    // generate another 50% more numbers
    int new_size_2 = unique_size * 1.5;
    NumericVector final_index_vector(new_size_2);

    // while there is not enough unique number
    while (unique_size < size){

      // generate the required amount of numbers
      NumericVector tmp_index_vector(new_size_2-unique_size);
      for(int i = 0 ; i < (new_size_2-unique_size) ; i++){
        tmp_index_vector[i] = unif_c(rng);
      }

      // put in "final_index_vector" the unique numbers (i.e. index_vector)
      // and the just generated additional numbers (i.e. tmp_index_vector) that may contain duplicated number
      std::copy(index_vector.begin(), index_vector.end(), final_index_vector.begin());
      std::copy(tmp_index_vector.begin(), tmp_index_vector.end(), final_index_vector.begin() + unique_size);

      // keep only unique numbers
      index_vector = unique(final_index_vector);
      // compute the amount of unique numbers
      unique_size = index_vector.length();
    }

    return(index_vector);
  }

}
