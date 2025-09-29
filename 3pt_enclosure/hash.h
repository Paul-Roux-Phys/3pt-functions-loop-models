#ifndef HASH_H_
#define HASH_H_

#include "common.h"

/* Make all entries of hash table point to NULL */
void init_hash()
{
  int i,j;
  
  for(i=0; i<hash_size; ++i)
    {
      for(j=0; j<2; ++j)
	p_hash[j][i] = NULL;
    }
}

/* Recursively empty a pointer list */
void eliminate(T_entry* p_entry)
{
  if(p_entry)
    {
      eliminate(p_entry->p_next);
      free(p_entry);
    }
}

/* Empty the entire hash table (at position data_pos) */
void empty_hash(int data_pos)
{
  int i;
  T_entry* p_entry;
  
  for(i=0; i<hash_size; ++i)
    {
      p_entry = p_hash[data_pos][i];
      eliminate(p_entry);
      p_hash[data_pos][i] = NULL;
    }
}

void construct_key()
{
  int i,k;
  int table[N+ORDINARY];
  unsigned char key0,flags;

  /* Order unmarked clusters */
  for(i=0; i<N+ORDINARY; i++)
    table[i]=0;
  k=ORDINARY-1;
  for(i=0; i<N; i++)
    {
      if(key[i]%32 >= ORDINARY)
	{
	  flags = key[i] & (7 << 5);
	  key0 = key[i]%32;
	  if(table[key0] != 0)
	    key[i]=table[key0] | flags;
	  else
	    {
	      k++;
	      table[key0] = k;
	      key[i] = k | flags;
	    }
	}
    }
}

void insert(int data_pos)
{
  int i,j;
  unsigned int hash;
  int is_weight_zero,no_strings;
  T_entry* p_entry;

  /* Don't insert something with zero weight */
  if(fabs(weight) < 1e-16)
    return;

  /* Count strings */
  no_strings = 0;
  for(i=0; i<N; ++i)
    if((key[i] == BOTTOM) || (key[i] == MIDDLE))
      ++no_strings;
  if(no_strings < min_strings)
    return;

  construct_key();
  
  // printf("inserting "); write(key);

  /* Compute hash */
  hash = 0;
  for (i=size-1; i>=0; --i)
    hash = (hash * hash_modulo + (unsigned char) key[i]) % hash_size;
  
  /* See if there is already an entry like that */
  p_entry = p_hash[data_pos][hash];
  while(p_entry)
    {
      if (Strcmp(p_entry->key,key)==0)
	{
	  /* Found! */
	  p_entry->weight += weight;
	  return;
	}
      p_entry = p_entry->p_next;
    }
  
  /* There wasn't. Thus insert a new entry */
  p_entry = (T_entry*)malloc(sizeof(T_entry));
  Strcpy(p_entry->key,key);
  p_entry->weight = weight;
  p_entry->p_next = p_hash[data_pos][hash];
  p_hash[data_pos][hash] = p_entry;
}

void insert_special(int data_pos)
{
  int i,j;
  unsigned int hash;
  int is_weight_zero;
  T_entry* p_entry;

  /* Don't insert something with zero weight */
  if(fabs(weight) < 1e-16)
    return;

  construct_key();
  
  //  printf("inserting "); write(key);

  /* Compute hash */
  hash = 0;
  for (i=size-1; i>=0; --i)
    hash = (hash * hash_modulo + (unsigned char) key[i]) % hash_size;
  
  /* See if there is already an entry like that */
  p_entry = p_hash[data_pos][hash];
  while(p_entry)
    {
      if (Strcmp(p_entry->key,key)==0)
	{
	  /* Found! */
	  p_entry->weight += weight;
	  return;
	}
      p_entry = p_entry->p_next;
    }
  
  /* There wasn't. Thus insert a new entry */
  p_entry = (T_entry*)malloc(sizeof(T_entry));
  Strcpy(p_entry->key,key);
  p_entry->weight = weight;
  p_entry->p_next = p_hash[data_pos][hash];
  p_hash[data_pos][hash] = p_entry;
}

double get_weight(int data_pos)
{
  int i,j;
  unsigned int hash;
  int is_weight_zero;
  T_entry* p_entry;

  construct_key();
  
  /* Compute hash */
  hash = 0;
  for (i=size-1; i>=0; --i)
    hash = (hash * hash_modulo + (unsigned char) key[i]) % hash_size;
  
  /* See if there is already an entry like that */
  p_entry = p_hash[data_pos][hash];
  while(p_entry)
    {
      if (Strcmp(p_entry->key,key)==0)
	{
	  /* Found! */
	  return p_entry->weight;
	}
      p_entry = p_entry->p_next;
    }
  
  /* There wasn't. */
  return 0.0;
}

void write_hash(int data_pos)
{
  int i,j;
  T_entry* p_entry;
  
  printf("HASH TABLE:\n");
  for(i=0; i<hash_size; ++i)
    {
      p_entry = p_hash[data_pos][i];
      while(p_entry)
	{
	  Strcpy(key,p_entry->key);
	  printf("weight=%.16f ",p_entry->weight);
	  printf("___ Key= ");
	  write(key);
	  p_entry = p_entry->p_next;
	}
    }
  printf("\n");
}

void write_hash_statistics(int data_pos)
{
  int i,j;
  int entries,min,max,level;
  T_entry* p_entry;
  
  entries = 0;
  min = -1;
  max = 0;
  
  for(i=0; i<hash_size; ++i)
    {
      p_entry = p_hash[data_pos][i];
      level = 0;
      while(p_entry)
	{
	  ++level;
	  p_entry = p_entry->p_next;
	}
      entries += level;
      if(min==-1) /* First time */
	min = level;
      if(level < min)
	min = level;
      if(level > max)
	max = level;
    }
  printf("Min, Max, Mean length: %d, %d, %d/%d = %.1f \n", 
	 min,max,entries,hash_size,(double)entries/hash_size);
  fflush(NULL);
}

void write(char *number)
{
  int i;
  
  printf("[ ");
  for(i=0; i<size; ++i)
    {
      printf("%d  ",number[i]);
    }
  printf("] \n");
}


void Strcpy(char *key1,char *key2)
{
  int i;
  
  for(i=0; i<size; ++i)
    key1[i] = key2[i];
}


int Strcmp(char *key1,char *key2)
{
  int i,result;
  
  result = 0;
  for(i=0; i<size; ++i)
    if(key1[i] != key2[i])
      return 1;
  return 0;
}

#endif // HASH_H_
