#ifndef ZEROBOARD_H
#define ZEROBOARD_H

typedef struct coord_set coord_set;
typedef struct table_item table_item;
typedef struct Board Board;

struct coord_set {
  int* coords;        // pointer to the indexes of the dataset values contained in this permutation of the query value
  int coords_size;    // number of indexes stored in coords (i.e. length of the permutation or no. of values making up the query value)
  coord_set* next;
};

struct table_item {
  // full value for the bin 
  int value;
  coord_set* head;
  coord_set* tail;
  // The variables prev and next are used to iterate through the table
  // contents in time: O(table bins used).
  // This 'linked-list' style of reaching numerically close combinations
  // is useful when catching all solutions within epsilon of the query value.
  int prev;
  int next;
  // This variable is used when there is a second table item that matches the bin
  // but has a different value. This can occur when there is a value difference
  // after the number of decimal places considered in the zeroboard, e.g.
  //    zb_dp = 3
  //    val_a = 162.10146432
  //    val_b = 162.10172009
  // It is used as a linked-list: simply a reference to the next table_item that
  // is stored in the same bin. It is a relatively rare occurance in cases where
  // the input set values are numerically distant.
  table_item* overflow;
};

struct Board{
  table_item* table;
  int table_size;     // calculated and then used for allocating memory and initialising start value and bin values 
  int start;          // start and end used for keeping track of the first and last used bins in the table
  int end;
};

Board* new_zeroboard(int digits, int dp_precision, int zb_dp) {
  // Size is relative to number of digits before the dp
  // and the number of values after the dp to successfully differentiate most values.
  // In this case, the number of digits before the dp is limited to 4 because masses
  // for Mass Spec won't go above the single thousands (ie. 9999) and the number of dp
  // accuracy these machines have is typically limited to 5, so we would expect that 
  // 3 dp is a good balance between distinguishing unique values in the zeroboard and
  // using the least memory possible. So that's 4 digits plus 3, which is 7, meaning
  // an array of 10million ints. This is merely one case however; the number of
  // decimal places accuracy can be defined by the user.
  Board* new_zb = malloc(sizeof(Board));
  new_zb->table_size = (int)pow(10,digits+zb_dp);
  new_zb->table = malloc(sizeof(table_item)*new_zb->table_size);
  new_zb->start = new_zb->table_size;
  new_zb->end = 0;
  for (int i=0; i<new_zb->table_size; ++i) {
    new_zb->table[i].value = -1;
    new_zb->table[i].next = -1;
  }
  return new_zb;
}

void board_insert(Board* zeroboard, double new_value, int* coords, int coords_size, int dp, int zb_dp) {
  // converting the query value to the representation used by the zeroboard storage table
  int bin_value = (new_value*pow(10,zb_dp));
  // if value not already inserted
  if (zeroboard->table[bin_value].value == -1) {
    zeroboard->table[bin_value].value = round(new_value*pow(10,dp));
    zeroboard->table[bin_value].overflow = NULL;

    coord_set* new_coords = malloc(sizeof(coord_set));
    new_coords->coords = coords;
    new_coords->coords_size = coords_size;
    new_coords->next = NULL;
    zeroboard->table[bin_value].head = new_coords;
    zeroboard->table[bin_value].tail = new_coords;
    // insert first item
    if (zeroboard->end == 0) {
      zeroboard->start = bin_value;
      zeroboard->end = bin_value;
      return;
    }
    // new start
    if (bin_value < zeroboard->start) {
      zeroboard->table[bin_value].next = zeroboard->start;
      zeroboard->table[zeroboard->start].prev = bin_value;
      zeroboard->start = bin_value;
      return;
    }
    // new end
    if (bin_value > zeroboard->end) {
      zeroboard->table[bin_value].prev = zeroboard->end;
      zeroboard->table[zeroboard->end].next = bin_value;
      zeroboard->end = bin_value;
      return;
    }
    // between start and end
    int counter = bin_value-1;
    while (zeroboard->table[counter].value == -1) {
      --counter;
    }
    zeroboard->table[bin_value].next = zeroboard->table[counter].next;
    zeroboard->table[bin_value].prev = counter;
    zeroboard->table[counter].next = bin_value;
    zeroboard->table[zeroboard->table[bin_value].next].prev = bin_value;
    return;
  
  // value already inserted (matches to user-defined number of dp)
  // so just add a set of coords to the list of coords for that value
  } else if ((zeroboard->table[bin_value].value) == (int)round(new_value*pow(10,dp))) {
    coord_set* new_coords = malloc(sizeof(coord_set));
    new_coords->coords = coords;
    new_coords->coords_size = coords_size;
    new_coords->next = NULL;
    // insert list item at head of list
    new_coords->next = zeroboard->table[bin_value].head;
    zeroboard->table[bin_value].head = new_coords;
    return;
  
  // value mismatch (to user-defined number of dp)
  // this must take into account the difference between the dp accuracy and the number of dp included in the table size
  // this must also be reflected in the get and contains methods
  } else {
    // write the new coord set for storing
    coord_set* new_coords   = malloc(sizeof(coord_set));
    new_coords->coords      = coords;
    new_coords->coords_size = coords_size;
    new_coords->next        = NULL;
    // check if value exists in the list
    table_item *item      = &zeroboard->table[bin_value],
               *last_item = item;
    while (item != NULL && (int)round(new_value*pow(10,dp)) != item->value) {
      item = item->overflow;
      if (item != NULL)
        last_item = item;
    }
    // value does not exist in overflow list, so add a new table item to the overflow list
    if (item == NULL) {
      last_item->overflow = malloc(sizeof(table_item));
      last_item = last_item->overflow;
      last_item->value = round(new_value*pow(10,dp));
      last_item->overflow = NULL;
      last_item->head = new_coords;
      last_item->tail = new_coords;
    // value does exist in the overflow list, so add the coord set to that item
    } else {
      new_coords->next = item->head;
      item->head = new_coords;
    }
  }
}

table_item* zeroboard_get_table_item(Board* zeroboard, double tareval, int zb_dp, int dp, double epsilon) {
  // negative values do not exist in the zeroboard
  if (tareval < 0.0)
    return NULL;
  // check the first table item this bin for a value match at the specified no. of dp and if needed 
  // advance through the overflow list to find the value
  int tareval_upp_bound = round((tareval+epsilon)*pow(10,dp)),
      tareval_low_bound = round((tareval-epsilon)*pow(10,dp)),
      bin_val           = tareval*pow(10,zb_dp),
      query_val         = round(tareval*pow(10,dp));
  if (zeroboard->table[bin_val].value != -1) {
    table_item* item = &zeroboard->table[bin_val];
    while (item != NULL && (item->value > tareval_upp_bound && item->value < tareval_low_bound))
      item = item->overflow;
    // if the item is null the value was not found in this bin (incl its overflow list)
    return item;
  } else {
    return NULL;
  }
}

void delete_zeroboard(Board* zeroboard) {
  // for each used bin in the table...
  for (int bin=zeroboard->start; bin!=-1; ) {
    // for each item associated with this bin...
    table_item* item = &zeroboard->table[bin];
    while (item != NULL) {
      // free all coord set memory associated with the item in this bin
      coord_set* temp = item->head;
      coord_set* prev = item->head;
      while (temp != NULL) {
        temp = temp->next;
        free(prev->coords);
        free(prev);
        prev = temp;
      }
      if (item != &zeroboard->table[bin]) {
        // free memory for all overflow items attached to this bin
        table_item* prev = item;
        item = item->overflow;
        free(prev);
      } else
        item = item->overflow;
    }
    // advance to next used bin in the table
    bin = zeroboard->table[bin].next;
  }
  // free memory for the table and zeroboard object
  free(zeroboard->table);
  free(zeroboard);
}

void print_zeroboard(Board* zeroboard) {
  // print each item in the zeroboard table
  FILE *f = fopen("zb.txt", "w");
  for (int i=zeroboard->start; i!=-1; ) {
    if (zeroboard->table[i].value != -1){
      fprintf(f, "%d: ", i);
      // print each item associated with bin
      table_item* item = &zeroboard->table[i];
      while (item != NULL) {
        fprintf(f, "%d -> ", item->value);
        // print all coord sets associated with each item
        coord_set* temp = item->head;
        while (temp != NULL) {
          for (int i=0; i<temp->coords_size; ++i)
            fprintf(f, "%d ", temp->coords[i]);
          fprintf(f, "; ");
          temp = temp->next;
        }
        item = item->overflow;
      }
      fprintf(f, "\n");
    }
    // advance to next item in the table
    i = zeroboard->table[i].next;
  }
  fclose(f);
}

#endif