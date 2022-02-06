#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include "zeroboard.h"

clock_t start, finish;
double cpu_time_used;
void quicksort();
void process_inputs();
void writeZeroBoard();
void queryZeroBoard();
void check_seq();
void get_results();
char known_mod;
char labels[] = {'G','A','S','P','V','T','C','I','N','D','Q','K','E','M','H','F','R','c','Y','W'};
char mods_lbl[] = {'0', '1', '2', '3', '4', '5'};

int main(int argc, char *argv[]) {

  // run-time variables
    double _dataset[] = {
      57.02146, 71.03711, 87.03203, 97.05276, 99.06841, 101.04768, 103.00919, 113.08406, 114.04293,
      115.02694, 128.05858, 128.09496, 129.04259, 131.04049, 137.05891, 147.06841, 156.10111, 160.03065, 163.06333, 186.07931
    };
    double mods[] = {
      -17.02655, -1.03164, 14.96327, 15.99491, 31.98982, 47.98473
    };

    int     dataset_size        = sizeof(_dataset)/sizeof(_dataset[0]),
            mods_size           = sizeof(mods)/sizeof(mods[0]),
            zb_size             = 13,    // sets a specific zeroboard size ; if 0 then automated
            zb_max              = 15,   // sets the maximum zeroboard size; if 0 then automated (a reasonable max is dependent on the size of the input dataset and number of queries to be performed)
            digits              = 4,    // max number of digits before the dp
            dp_precision        = 4,    // total number of dp precision
            zeroboard_precision = 4,    // number of dp considered in the zeroboard
            combination_length  = 0;    // if specified, algorithm checks that combination length only
    double *dataset             = &_dataset[0],
            epsilon             = 0.0002;  // the amount by which the query value can vary

  // read query value from command line argument and error check it
    double query_value;
    sscanf(argv[1], "%lf", &query_value);
    sscanf(argv[2], "%c", &known_mod);
  
  // error check input values, sort input dataset, and remove duplicates
    process_inputs(query_value, zb_size, zb_max, digits, dp_precision, zeroboard_precision, combination_length, epsilon, dataset, &dataset_size);
    
    if (query_value < dataset[0]) {
      printf("\nERROR: Query value cannot be less than dataset minimum\n\tQuery value: %f\n\tDataset min: %f\n\n", query_value, dataset[0]);
      exit(EXIT_FAILURE);
    }
    printf("Query Value: %.5f\n", query_value);

  // create the zeroboard
    start              = clock();
      Board* zeroboard = new_zeroboard(digits, dp_precision, zeroboard_precision);
      writeZeroBoard(dataset, dataset_size, zb_size, zb_max, dp_precision, zeroboard_precision, query_value, zeroboard);
      //print_zeroboard(zeroboard);
    finish             = clock();
    cpu_time_used      = ((double) (finish - start)) / CLOCKS_PER_SEC;
    printf("\n%f seconds to create zeroboard\n\n", cpu_time_used);

  // query the zeroboard
    start         = clock();
      queryZeroBoard(dataset, dataset_size, &mods[0], mods_size, zb_size, zb_max, dp_precision, zeroboard_precision, query_value, zeroboard, combination_length, epsilon);
    finish        = clock();
    cpu_time_used = ((double) (finish - start)) / CLOCKS_PER_SEC;
    printf("%f seconds to query zeroboard\n\n", cpu_time_used);

  // free heap memory used by the zeroboard
    start         = clock();
    delete_zeroboard(zeroboard);
    finish        = clock();
    cpu_time_used = ((double) (finish - start)) / CLOCKS_PER_SEC;
    printf("\n%f seconds to free zeroboard memory\n", cpu_time_used);
  
  return 0;
}

void process_inputs(double query_val, int zb_size, int zb_max, int digits, int dp, int zb_dp, 
                    int combin_len, double epsilon, double *S, int *size_S) {
  // PART 1: error check input values
    if (digits < 1 || dp < 0 || zb_dp < 0 || combin_len < 0 || epsilon < 0.0) {
      printf("\nERROR: Input values cannot be negative and number of digits must be at least 1\n"
      "\tDigits\t\t: %d\n"
      "\tDP\t\t: %d\n"
      "\tZeroboard DP\t: %d\n"
      "\tCombination Len\t: %d \t\t(if 0, all valid lengths are checked)\n"
      "\tEpsilon\t\t: %f \t(if 0, only precise solutions are included)\n\n",
      digits, dp, zb_dp, combin_len, epsilon);
      exit(EXIT_FAILURE);
    }
    if (zb_dp > dp) {
      printf("ERROR: Zeroboard precision cannot be greater than decimal places precision\n");
      printf("\tzeroboard precision: %d\n\tdp precision: %d\n\n", zb_dp, dp);
      exit(EXIT_FAILURE);
    }
    if (zb_size < 0 || zb_max < 0) {
      printf("ERROR: Zeroboard size and zeroboard max must be >= 0\n");
      printf("\tzeroboard size: %d\n\tzeroboard max : %d\n\n", zb_size, zb_max);
      exit(EXIT_FAILURE);
    }

  // PART 2: If input is not sorted, sort it
    for (int i=1; i<*size_S; ++i)
      if (S[i-1] > S[i])
        quicksort(S, 0, *size_S);

  // PART 3: use sorted characteristic to remove duplicates in linear time
    bool duplicates = false;
    // check for duplicates
    for (int i=0; i<(*size_S)-1; ++i) {
      // if value is duplicated, set variable to true and break from iteration
      if (S[i] == S[i+1]) {
        duplicates == true;
        break;
      }
    }
    // if duplicates exist, remove them
    if (duplicates) {
      // create new dataset
      double *remove_duplicates = malloc(sizeof(double)*(*size_S));
      int counter = 0;
      // copy unqiue values to new dataset 
      for (int i=0; i<*size_S; ) {
        remove_duplicates[i] = S[i];
        counter = i+1;
        while (S[counter] == S[i])
          ++counter;
        i = counter;
      }
      // redirect pointer and reset dataset size
      free(S);
      S = remove_duplicates;
      *size_S = counter;
    }

}

void writeZeroBoard(double *S, int n, int zb_size, int zb_max, int dp, int zb_dp, double query_val, Board* zeroboard) {
  // method variables: size of zeroboard and size of tracking array
    int n_zerobased = n -1,
        min_zb_size     = 3,
        min_len         = (int)(query_val/S[n-1]);
    if (min_len < min_zb_size)
        min_len         = min_zb_size;
    if (zb_max != 0 && min_len > zb_max)
        min_len         = zb_max;
    if (zb_size != 0)
        min_len         = zb_size;
    //printf("zb size: %d\n", min_len);
    int zb_tracker_len  = min_len - min_zb_size + 1,
        zb_tracker[zb_tracker_len];
    // initialise zeroboard tracker array values
    for (int i=0; i<zb_tracker_len; ++i)
        zb_tracker[i]   = 0;
    double  dec_places  = pow(10,dp),
            S_max       = (S[n_zerobased]);
    printf("Zeroboard size: %d\n", min_len);
       
  // iterate through tracking array (which maintains tracking of current zeroboard triangle)
    while (zb_tracker[zb_tracker_len-1] <= n_zerobased) {
      while (zb_tracker[0] <= n_zerobased) {
        // iterate through columns and rows of current zeroboard triangle
        int colCounter = zb_tracker[0];
        while (colCounter <= n_zerobased) {
          int rowCounter = colCounter;
          while (rowCounter <= n_zerobased) {
            // calculate new value for insertion according to position in zeroboard (by tracker array and row/col counters)
            double new_value = (S_max - ((S[colCounter]))) + (S_max - (S[rowCounter]));
            int coord_list_size = zb_tracker_len + 2,
                counter         = 0,
                *coord_list     = malloc(sizeof(int)*coord_list_size);
            coord_list[coord_list_size-1] = rowCounter;
            coord_list[coord_list_size-2] = colCounter;
            while (counter < zb_tracker_len) {
              new_value = new_value + (S_max - (S[zb_tracker[counter]]));
              coord_list[coord_list_size - counter - 3] = zb_tracker[counter];
              ++counter;
            }
            // calculate tare value (numberic distance from zeroboard max) and insert into zeroboard
            board_insert(zeroboard, new_value, coord_list, coord_list_size, dp, zb_dp);
            ++rowCounter;
          }
          ++colCounter;
        }
        ++zb_tracker[0];
      }

      // maintain tracking array
      int decrementor = 0;
      while (zb_tracker[decrementor] >= n-1 && decrementor < zb_tracker_len-1)
        ++decrementor;
      ++zb_tracker[decrementor];
      while (decrementor > 0) {
        --decrementor;
        zb_tracker[decrementor] = zb_tracker[decrementor+1];
      }
    }
}

void queryZeroBoard(bool *soln_found, char *known_seq, double *S, int n, double *mods, int mods_n, int _zb_size, int zb_max, int dp, int zb_dp, double query_val, Board* zeroboard, int combination_length, double epsilon) {
  // method variables
    // min zeroboard size hard coded to permutation length 3
    int     min_zb_size     = 3,
            zb_size         = (int)(query_val/S[n-1]);
    // if the zb size is less than the min, set it to the min which is 3
    if (zb_size < min_zb_size)
            zb_size         = min_zb_size;
    // if the zb max size is specified by the user and the calculated size is greater, set it to the user specified value
    if (zb_max != 0 && zb_size > zb_max)
            zb_size         = zb_max;
    // if the user specifies the zb size, set it to that size
    int     end_length      = zb_size;
    if (_zb_size != 0) {
            zb_size         = _zb_size;
            end_length      = _zb_size;
    }
    int     acidCounter     = (int)(query_val/S[0]), 
            n_zeroBased     = n-1;
    ssize_t resultsCounter  = 0,
            totalResults    = 0;
    double  tareMass        = 0.0,
            dec_places      = pow(10,dp),
            S_max           = (S[n_zeroBased]),
            query_val_dp    = query_val;

    // if combination length set, only search that length
    if (combination_length != 0) {
      acidCounter           = combination_length;
      end_length            = combination_length-1;
    }
    
    // tracking arrays
    int     array_size      = acidCounter-zb_size,
            array[array_size];
    double  mins[array_size],
            maxs[array_size],
    // set key sequences
            keySeq[n],
            keySeqSum[n];
    for (int i=0; i<n; ++i) {
      if (i == 0) {
        keySeq[i]     = 0.0;
        keySeqSum[i]  = 0.0;
      } else {
        keySeq[i]     = S[i] - S[i-1];
        keySeqSum[i]  = keySeqSum[i-1] + keySeq[i];
      }
    }

  printf("combination length : results %d\n", end_length);
  // iterate through valid combination lengths above minimum
  while (acidCounter > end_length && acidCounter*S[n_zeroBased] >= query_val) {
    int     dim     = 0;
    double  acidMax = acidCounter*S[n_zeroBased],
            acidMin = acidCounter*S[0],
         acidMax_dp = acidMax;
    
    if ((int)acidMax_dp == (int)query_val_dp) {
      // the coord set will be the set max to the current combination length
      ++resultsCounter;
    } else if ((int)(acidMin) == (int)query_val_dp) {
      // the coord set will be the set min to the current combination length
      ++resultsCounter;
    } else {
      // set tracking arrays
        for (int i=0; i<array_size; ++i) {
          mins[i] = S[0]*acidCounter;
          maxs[i] = S[0]*(i+1) + S[n-1]*(acidCounter-(i+1));
          array[i] = 0;
        }
      while (mins[0] <= query_val && dim < acidCounter-zb_size) {
        // max finding
          while (dim <= acidCounter-zb_size && array[dim] < n_zeroBased) {
            while (maxs[dim] <= query_val && array[dim] < n_zeroBased) {
              // check zeroboard if the maxs value is the query value
              if ((int)maxs[dim]*(int)dec_places == (int)query_val*(int)dec_places) {
                int u = 0;
                double tareMass = 0.0;
                while (u < acidCounter - zb_size) {
                  tareMass = tareMass - ((S_max)-(S[array[u]]));
                  ++u;
                }
                tareMass = (tareMass + ((acidMax_dp - query_val_dp)));
                // check using calculated tare mass
                get_results(S, mods, mods_n, zeroboard, zb_size, query_val_dp, tareMass, &resultsCounter, array, array_size, acidCounter-zb_size-1, zb_dp, dp, epsilon);
              }
              ++array[dim];
              maxs[dim] = maxs[dim] + (S[array[dim]]-S[array[dim]-1]);
            }

            if (dim < acidCounter - zb_size - 1 && array[dim] < n_zeroBased) {
              ++dim;
              array[dim] = array[dim-1];
              maxs[dim] = maxs[dim-1] - (S[n_zeroBased]-S[array[dim-1]]);
            } else
              break;
          }
        // min finding
          if (array[dim] <= n_zeroBased) {
            mins[dim] = maxs[dim] - (S[n_zeroBased]-S[array[dim]])*zb_size;
            while (mins[dim] <= query_val && array[dim] < n_zeroBased) {
              // check zeroboard
              int u = 0;
              double tareMass = 0.0;
              while (u < acidCounter - zb_size) {
                tareMass = tareMass - (S_max-S[array[u]]);
                ++u;
              }
              tareMass = ((tareMass + (acidMax_dp - query_val_dp)));
              get_results(S, mods, mods_n, zeroboard, zb_size, query_val_dp, tareMass, &resultsCounter, array, array_size, acidCounter-zb_size-1, zb_dp, dp, epsilon);
              ++array[dim];
              mins[dim] = mins[dim] + (S[array[dim]]-S[array[dim]-1])*(zb_size+1);
            }
          }
        // maintain tracker array
          while (mins[dim] > query_val && dim > 0) {
            if (dim > 0) {
              --dim;
              ++array[dim];
              double temp = 0;
              int u = 0;
              while (u < dim) {
                temp = temp + S[array[u]];
                ++u;
              }
              mins[dim] = temp + S[array[dim]]*(acidCounter-dim);
            } else if (dim == 0)
              mins[dim] = acidMin + keySeqSum[array[dim]]*acidCounter;
          }
        // reset tracker array lower dimensions
          int r = dim+1;
          while (r < acidCounter-zb_size) {
            array[r] = array[r-1];
            ++r;
          }
          ++dim;
          maxs[dim] = mins[dim-1] + (S[n_zeroBased]-S[array[dim]]) * (acidCounter-(dim+1));
      }
    }
    
    printf("\t%d\t\t%ld\n", acidCounter, resultsCounter);
    totalResults = totalResults + resultsCounter;
    resultsCounter = 0;
    --acidCounter;
  }

  // check min combination length (the zeroboard size)
  if (combination_length == 0) {
    tareMass = (acidCounter*S[n_zeroBased] - query_val);
    get_results(soln_found, known_seq, S, zeroboard, zb_size, query_val_dp, tareMass, &resultsCounter, array, array_size, acidCounter-zb_size-1, zb_dp, dp, epsilon);
    printf("\t%d\t\t%ld\n", acidCounter, resultsCounter);
    totalResults = totalResults + resultsCounter;
  }

  printf("results: %ld\n", totalResults);
}

void get_results(double *dataset, double *mods, int mods_n, Board* zeroboard, int zb_size, double query_val_dp, double taremass, int* num_results, int* array, int array_size, int combin_len, int zb_dp, int dp, double epsilon) {
  table_item* ptr = zeroboard_get_table_item(zeroboard, taremass, zb_dp, dp, epsilon);
  // function variables
    int          n = combin_len+zb_size+1;
    coord_set   *set;
    char        *this_seq = malloc(sizeof(char)*(n));
  // All combinations for precise solution included - gather all combinations in list for this value
    if (ptr != NULL) {
      while (ptr != NULL) {
        set = ptr->head;
        while (set != NULL) {
          if (set->coords[0] >= array[combin_len] || combin_len == -1) {
            // print sequence
            double val = 0.0;
            for (int i=0; i<=combin_len; ++i) {
              val = val + dataset[array[i]];
              this_seq[i] = labels[array[i]];
            }
            for (int i=0; i<set->coords_size; ++i) {
              val = val + dataset[set->coords[i]];
              this_seq[i+combin_len+1] = labels[set->coords[i]];
            }
            int val_int = round(val*pow(10,dp));
            if (val_int >= (int)round((query_val_dp-epsilon)*pow(10,dp)) && val_int <= (int)round((query_val_dp+epsilon)*pow(10,dp))) {
              for (int i=0; i<n; ++i)
                printf("%c", labels[this_seq[i]]);
              printf("\n");
              ++(*num_results);
            }
          } else
            break;
          set = set->next;
        }
        ptr = ptr->overflow;
      }
    }

    // If there is a modification:
      if (known_mod != 'x')
      // For every modification...
        for (int i=0; i<mods_n; ++i) {
          // ...calculate the taremass and check if that is in the zeroboard
          int bin = (int)((taremass+mods[i])*pow(10,zb_dp));
          if (bin >=0 && bin < zeroboard->table_size && zeroboard->table[bin].value != -1) {
            ptr = &zeroboard->table[bin];
            while (ptr != NULL) {
              set = ptr->head;
              while (set != NULL) {
                if (set->coords[0] >= array[combin_len] || combin_len == -1) {
                  // record coords
                  double val = mods[i];
                  for (int i=0; i<=combin_len; ++i) {
                    val = val + dataset[array[i]];
                    this_seq[i] = labels[array[i]];
                  }
                  for (int i=0; i<set->coords_size; ++i) {
                    val = val + dataset[set->coords[i]];
                    this_seq[i+combin_len+1] = labels[set->coords[i]];
                  }
                  int val_int = (val*pow(10,dp));
                  if (val_int >= (int)round((query_val_dp-epsilon)*pow(10,dp)) && val_int <= (int)round((query_val_dp+epsilon)*pow(10,dp))) {
                    printf("%c", mods_lbl[i]);
                    for (int i=0; i<n; ++i)
                      printf("%c", labels[this_seq[i]]);
                    printf("\n");
                    ++(*num_results);
                  }
                } else
                  break;
                set = set->next;
              }
              ptr = ptr->overflow;
            }
          }
        }

  // if epsilon is specified, get all values within range of (query_val-epsilon) and (query_val+epsilon)
    if (epsilon != 0.0) {
    // calculate bounds
      int upper_bound = round((taremass+epsilon)*pow(10,dp)),
          lower_bound = ((taremass-epsilon)*pow(10,dp)),
          upper_bin   = round((taremass+epsilon)*pow(10,zb_dp)),
          lower_bin   = ((taremass-epsilon)*pow(10,zb_dp)),
          bin_val     = (taremass*pow(10,zb_dp))+1;
    // if the pointer is null, the bin was empty so advance to the next used bin that is within epsilon range of the query value...
      if (bin_val > 0 && bin_val < zeroboard->table_size) {
      // get next bucket
        table_item* upper_range;
        while (zeroboard->table[bin_val].value == -1 && bin_val <= upper_bin)
          ++bin_val;
        upper_range = &zeroboard->table[bin_val];
      // if within range, move forward in zeroboard storage table gathering combinations with sum less than (query_val+epsilon)
        if (upper_range->value != -1) {
          while (upper_range->value <= upper_bound) {
            // check all items in this bin, including overflow
            table_item *bin_item = upper_range;
            while (bin_item != NULL) {
              if (bin_item->value <= upper_bound) {
                // gather all combinations in list for this value
                set = bin_item->head;
                while (set != NULL) {
                  if (set->coords[0] >= array[combin_len] || combin_len == -1) {
                    double val = 0.0;
                    for (int i=0; i<=combin_len; ++i) {
                      val = val + dataset[array[i]];
                      this_seq[i] = labels[array[i]];
                    }
                    for (int i=0; i<set->coords_size; ++i) {
                      val = val + dataset[set->coords[i]];
                      this_seq[i+combin_len+1] = labels[set->coords[i]];
                    }
                    if ((int)round(val*pow(10,dp)) <= (int)round((query_val_dp+epsilon)*pow(10,dp))) {
                      for (int i=0; i<n; ++i)
                          printf("%c", labels[this_seq[i]]);
                      printf("\n");
                      ++(*num_results);
                    }
                  } else
                    break;
                  set = set->next;
                }
              }
              bin_item = bin_item->overflow;
            }
            // if next item exists, move pointer to address of next item
            if (upper_range->next != -1)
              upper_range = &zeroboard->table[upper_range->next];
            else
              break;
          }
        }
      }

    // move to previous bucket in zeroboard storage table with value greater than query value minus epsilon
      // get previous bucket
        table_item* lower_range;
        bin_val = taremass*pow(10,zb_dp)-1;
        if (bin_val > 0) {
          while (zeroboard->table[bin_val].value == -1 && bin_val >= lower_bin)
            --bin_val;
          lower_range = &zeroboard->table[bin_val];
          // move backwards in zeroboard storage table gathering combinations with sum greater than (query_val-epsilon)
          if (lower_range->value != -1) {
            while (lower_range->value >= lower_bound) {
              // check all items in this bin, including overflow
              table_item *bin_item = lower_range;
              while (bin_item != NULL) {
                if (bin_item->value >= lower_bound) {
                  // gather all combinations in list for this value
                  set = bin_item->head;
                  while (set != NULL) {
                    if (set->coords[0] >= array[combin_len] || combin_len == -1) {
                      double val = 0.0;
                      for (int i=0; i<=combin_len; ++i) {
                        val = val + dataset[array[i]];
                        this_seq[i] = labels[array[i]];
                      }
                      for (int i=0; i<set->coords_size; ++i) {
                        val = val + dataset[set->coords[i]];
                        this_seq[i+combin_len+1] = labels[set->coords[i]];
                      }
                      if ((int)round(val*pow(10,dp)) >= (int)round((query_val_dp-epsilon)*pow(10,dp))) {
                        for (int i=0; i<n; ++i)
                          printf("%c", labels[this_seq[i]]);
                        printf("\n");
                        ++(*num_results);
                      }
                    } else
                      break;
                    set = set->next;
                  }
                }
                bin_item = bin_item->overflow;
              }
              // if previous item exists, move pointer to address of previous item
              if (lower_range->prev != -1)
                lower_range = &zeroboard->table[lower_range->prev];
              else
                break;
            }
          }
        }
    }
  free(this_seq);
}

void quicksort(double a[], int first, int last) {
  int pivot, i, j;
  if(first < last) {
    pivot = first;
    i = first;
    j = last;
    while (i < j) {
      while(a[i] <= a[pivot] && i < last)
        i++;
      while(a[j] > a[pivot])
        j--;
      if(i < j) {
        double temp;
        temp = a[i];
        a[i] = a[j];
        a[j] = temp;
      }
    }
    double temp;
    temp = a[pivot];
    a[pivot] = a[j];
    a[j] = temp;
    quicksort(a, first, j - 1);
    quicksort(a, j + 1, last);
  }
}
