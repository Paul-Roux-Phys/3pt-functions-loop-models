#ifndef TRANSFER_H_
#define TRANSFER_H_
#include "common.h"

double vector_norm(int data_pos)
{
  int i;
  double square_sum;
  T_entry* p_entry;
  
  square_sum = 0.0;
  for(i=0; i<hash_size; ++i)
    {
      p_entry = p_hash[data_pos][i];
      while(p_entry)
	{
	  square_sum += (p_entry->weight) * (p_entry->weight);
	  p_entry = p_entry->p_next;
	}
    }
  return sqrt(square_sum);
}

void normalise_vector(int data_pos)
{
  int i;
  double v_norm;
  T_entry* p_entry;

  v_norm = vector_norm(data_pos);
  exponent += log(v_norm);
  for(i=0; i<hash_size; ++i)
    {
      p_entry = p_hash[data_pos][i];
      while(p_entry)
	{
	  p_entry->weight /= v_norm;
	  p_entry = p_entry->p_next;
	}
    }
}

void unpack_state(T_entry* p_entry)
{
  int i;
  
  weight = p_entry->weight;
  Strcpy(key,p_entry->key);
}

void join(int pos1,int pos2)
{
  int p1,p2;
  int flag1,flag2,flag,val;

  flag1 = (key[pos1] >> 5) & 7;
  flag2 = (key[pos2] >> 5) & 7;
  flag = flag1 ^ flag2;
  if((key[pos1]%32 >= ORDINARY) && (key[pos2]%32 >= ORDINARY))
    {
      for(p1=0; (key[p1] != key[pos1]) || (p1 == pos1); ++p1);
      for(p2=0; (key[p2] != key[pos2]) || (p2 == pos2); ++p2);
      val = key[p1]%32;
      key[p1] = val | (flag << 5);
      key[p2] = val | (flag << 5);
    }
  else if(key[pos1]%32 >= ORDINARY)
    {
      for(p1=0; (key[p1] != key[pos1]) || (p1 == pos1); ++p1);
      val = key[pos2]%32;
      key[p1] = val | (flag << 5);
    }
  else if(key[pos2]%32 >= ORDINARY)
    {
      for(p2=0; (key[p2] != key[pos2]) || (p2 == pos2); ++p2);
      val = key[pos1]%32;
      key[p2] = val | (flag << 5);
    }
}

double getloopweight(unsigned char k)
{
  int flag;

  flag = (k >> 5) & 7;
  if(flag == 0)
    return n_loop;
  else if(flag == 1)
    return n2;
  else
    exit(1);
}


void toggle_seam(int pos,int seamno)
{
  unsigned char mask;
  int p;

  mask = 1 << (5 + seamno);
  key[pos] ^= mask;
  if(key[pos]%32 >= ORDINARY)
    {
      for(p=0; (key[p]%32 != key[pos]%32) || (p == pos); ++p);
      key[p] ^= mask;
    }
}

void insert_aux_space()
{
  int hash_entry,i;
  T_entry* p_entry;
  T_entry* temp;

#ifdef DEBUG
  printf("before insert_aux_space : ");
  write_hash(data_position);
#endif

  for(hash_entry=0; hash_entry<hash_size; ++hash_entry)
    {
      p_entry = p_hash[data_position][hash_entry];
      while(p_entry)
        {
          unpack_state(p_entry);

          /* Shift indices */
          for(i=n; i>0; --i)
            key[i] = key[i-1];

          key[0] = EMPTY;
          key[n+1] = EMPTY;
          insert(1-data_position);

          key[0] = ORDINARY+n+1;
          key[n+1] = ORDINARY+n+1;
          insert(1-data_position);
          
          temp=p_entry->p_next;
          free(p_entry);
          p_entry = temp;
        }
      p_hash[data_position][hash_entry]=NULL;
    }
  data_position = 1-data_position;

#ifdef DEBUG
  printf("after insert_aux_space : ");
  write_hash(data_position);
#endif
}

void insert_aux_space_special(int seam)
{
  int hash_entry,i;
  T_entry* p_entry;
  T_entry* temp;

#ifdef DEBUG
  printf("before insert_aux_space_special : ");
  write_hash(data_position);
#endif

  for(hash_entry=0; hash_entry<hash_size; ++hash_entry)
    {
      p_entry = p_hash[data_position][hash_entry];
      while(p_entry)
        {
          unpack_state(p_entry);

          /* Shift indices */
          for(i=n; i>0; --i)
            key[i] = key[i-1];

          key[0] = EMPTY;
          key[n+1] = EMPTY;
          insert_special(1-data_position);
	  
          key[0] = ORDINARY+n+1;
          key[n+1] = ORDINARY+n+1;
	  if(seam == 1)
	    toggle_seam(0,0);
          insert_special(1-data_position);
	  
          temp=p_entry->p_next;
          free(p_entry);
          p_entry = temp;
        }
      p_hash[data_position][hash_entry]=NULL;
    }
  data_position = 1-data_position;

#ifdef DEBUG
  printf("after insert_aux_space_special : ");
  write_hash(data_position);
#endif
}

void R_matrix(int pos)
{
  int hash_entry,i;
  int pos1;
  double save_weight;
  T_entry* p_entry;
  T_entry* temp;

#ifdef DEBUG
  printf("before R_matrix(%d) : ",pos);
  write_hash(data_position);
#endif

  pos1 = pos+1;
  
  for(hash_entry=0; hash_entry<hash_size; ++hash_entry)
    {
      p_entry = p_hash[data_position][hash_entry];
      while(p_entry)
        {
          unpack_state(p_entry);
	  save_weight = weight;
	  
          if((key[pos] == EMPTY) && (key[pos1] == EMPTY))
            {
              /* rho(1) delta(a,b,c,d) */
	      weight = save_weight * rho[1];
              insert(1-data_position);

              /* rho(5) delta(a,b,d) adj(a,c) */
              key[pos] = ORDINARY+n+1;
              key[pos1] = ORDINARY+n+1;
	      weight = save_weight * rho[5];
              insert(1-data_position);
            }

          else if((key[pos] == EMPTY) && (key[pos1] != EMPTY))
            {
              /* rho(3) delta(a,c,d) adj(a,b) */
	      weight = save_weight * rho[3];
              insert(1-data_position);
              
              /* rho(7) delta(a,d) delta(b,c) adj(a,b) */
              key[pos] = key[pos1];
              key[pos1] = EMPTY;
	      weight = save_weight * rho[7];
              insert(1-data_position);
            }

          else if((key[pos] != EMPTY) && (key[pos1] == EMPTY))
            {
              /* rho(2) delta(a,b,c) adj(a,d) */
	      weight = save_weight * rho[2];
              insert(1-data_position);
              
              /* rho(6) delta(a,b) delta(c,d) adj(a,c) */
              key[pos1] = key[pos];
              key[pos] = EMPTY;
	      weight = save_weight * rho[6];
              insert(1-data_position);
            }

          else if((key[pos] != EMPTY) && (key[pos1] != EMPTY))
            {
              /* rho(8) delta(a,c) adj(a,b) adj(a,d) */
	      weight = save_weight * rho[8];
              insert(1-data_position);

	      /* TL generator */
	      if((key[pos] >= ORDINARY) || (key[pos1] >= ORDINARY))
		{
		  if(key[pos] == key[pos1])
		    {
		      /* rho(4) delta(b,c,d) adj(a,b) */
		      key[pos] = EMPTY;
		      key[pos1] = EMPTY;
		      weight = save_weight * rho[4] * n_loop;
		      insert(1-data_position);
                  
		      /* rho(9) delta(b,d) adj(b,c) adj(a,b) */
		      key[pos] = ORDINARY+n+1;
		      key[pos1] = ORDINARY+n+1;
		      weight = save_weight * rho[9] * n_loop;
		      insert(1-data_position);
		    }
		  else
		    {
		      join(pos,pos1);

		      /* rho(4) delta(b,c,d) adj(a,b) */
		      key[pos] = EMPTY;
		      key[pos1] = EMPTY;
		      weight = save_weight * rho[4];
		      insert(1-data_position);
		      
		      /* rho(9) delta(b,d) adj(b,c) adj(a,b) */
		      key[pos] = ORDINARY+n+1;
		      key[pos1] = ORDINARY+n+1;
		      weight = save_weight * rho[9];
		      insert(1-data_position);
		    }
		}
	      else if(key[pos] != key[pos1])
		{
		  /* rho(4) delta(b,c,d) adj(a,b) */
		  key[pos] = EMPTY;
		  key[pos1] = EMPTY;
		  weight = save_weight * rho[4];
		  insert(1-data_position);
		      
		  /* rho(9) delta(b,d) adj(b,c) adj(a,b) */
		  key[pos] = ORDINARY+n+1;
		  key[pos1] = ORDINARY+n+1;
		  weight = save_weight * rho[9];
		  insert(1-data_position);
		}
            }
	  
          temp=p_entry->p_next;
          free(p_entry);
          p_entry = temp;
        }
      p_hash[data_position][hash_entry]=NULL;
    }
  data_position = 1-data_position;
  
#ifdef DEBUG
  printf("after R_matrix(%d) : ",pos);
  write_hash(data_position);
#endif
}

void R_matrix_special(int pos)
{
  int hash_entry,i;
  int pos1;
  double save_weight,loopfactor;
  T_entry* p_entry;
  T_entry* temp;

#ifdef DEBUG
  printf("before R_matrix_special(%d) : ",pos);
  write_hash(data_position);
#endif

  pos1 = pos+1;
  
  for(hash_entry=0; hash_entry<hash_size; ++hash_entry)
    {
      p_entry = p_hash[data_position][hash_entry];
      while(p_entry)
        {
          unpack_state(p_entry);
	  save_weight = weight;
	  // printf("treating "); write(key);
	  
          if((key[pos] == EMPTY) && (key[pos1] == EMPTY))
            {
              /* rho(1) delta(a,b,c,d) */
	      weight = save_weight * rho[1];
              insert_special(1-data_position);

              /* rho(5) delta(a,b,d) adj(a,c) */
              key[pos] = ORDINARY+n+1;
              key[pos1] = ORDINARY+n+1;
	      weight = save_weight * rho[5];
              insert_special(1-data_position);
            }

          else if((key[pos] == EMPTY) && (key[pos1] != EMPTY))
            {
              /* rho(3) delta(a,c,d) adj(a,b) */
	      weight = save_weight * rho[3];
              insert_special(1-data_position);
              
              /* rho(7) delta(a,d) delta(b,c) adj(a,b) */
              key[pos] = key[pos1];
              key[pos1] = EMPTY;
	      weight = save_weight * rho[7];
              insert_special(1-data_position);
            }

          else if((key[pos] != EMPTY) && (key[pos1] == EMPTY))
            {
              /* rho(2) delta(a,b,c) adj(a,d) */
	      weight = save_weight * rho[2];
              insert_special(1-data_position);
              
              /* rho(6) delta(a,b) delta(c,d) adj(a,c) */
              key[pos1] = key[pos];
              key[pos] = EMPTY;
	      weight = save_weight * rho[6];
              insert_special(1-data_position);
            }

         else if((key[pos] != EMPTY) && (key[pos1] != EMPTY))
            {
              /* rho(8) delta(a,c) adj(a,b) adj(a,d) */
	      weight = save_weight * rho[8];
              insert_special(1-data_position);

	      /* TL generator */
	      if((key[pos]%32 >= ORDINARY) || (key[pos1]%32 >= ORDINARY))
		{
		  if(key[pos]%32 == key[pos1]%32)
		    {
		      loopfactor = getloopweight(key[pos]);
		      
		      /* rho(4) delta(b,c,d) adj(a,b) */
		      key[pos] = EMPTY;
		      key[pos1] = EMPTY;
		      weight = save_weight * rho[4] * loopfactor;
		      insert_special(1-data_position);
                  
		      /* rho(9) delta(b,d) adj(b,c) adj(a,b) */
		      key[pos] = ORDINARY+n+1;
		      key[pos1] = ORDINARY+n+1;
		      weight = save_weight * rho[9] * loopfactor;
		      insert_special(1-data_position);
		    }
		  else
		    {
		      join(pos,pos1);
		      
		      /* rho(4) delta(b,c,d) adj(a,b) */
		      key[pos] = EMPTY;
		      key[pos1] = EMPTY;
		      weight = save_weight * rho[4];
		      insert_special(1-data_position);
		      
		      /* rho(9) delta(b,d) adj(b,c) adj(a,b) */
		      key[pos] = ORDINARY+n+1;
		      key[pos1] = ORDINARY+n+1;
		      weight = save_weight * rho[9];
		      insert_special(1-data_position);
		    }
		}
	      /* Connect strings if they have matching target codes and stay clear of seam */
	      else if(((target[key[pos]%32]==key[pos1]%32) && (target[key[pos1]%32]==key[pos]%32)) &&
		      ((key[pos] >> 5) == (key[pos1] >> 5)))
		{
		  /* rho(4) delta(b,c,d) adj(a,b) */
		  key[pos] = EMPTY;
		  key[pos1] = EMPTY;
		  weight = save_weight * rho[4];
		  insert_special(1-data_position);
		      
		  /* rho(9) delta(b,d) adj(b,c) adj(a,b) */
		  key[pos] = ORDINARY+n+1;
		  key[pos1] = ORDINARY+n+1;
		  weight = save_weight * rho[9];
		  insert_special(1-data_position);
		}
            }
	  
          temp=p_entry->p_next;
          free(p_entry);
          p_entry = temp;
        }
      p_hash[data_position][hash_entry]=NULL;
    }
  data_position = 1-data_position;
  
#ifdef DEBUG
  printf("after R_matrix_special(%d) : ",pos);
  write_hash(data_position);
#endif
}

void remove_aux_space()
{
  int hash_entry,i;
  int pos,pos1;
  T_entry* p_entry;
  T_entry* temp;

#ifdef DEBUG
  printf("before remove_aux_space : ");
  write_hash(data_position);
#endif

  pos = n;
  pos1 = n+1;
  
  for(hash_entry=0; hash_entry<hash_size; ++hash_entry)
    {
      p_entry = p_hash[data_position][hash_entry];
      while(p_entry)
        {
          unpack_state(p_entry);

          if((key[pos] == EMPTY) && (key[pos1] == EMPTY))
            insert(1-data_position);
          else if((key[pos] != EMPTY) && (key[pos1] != EMPTY))
	    if((key[pos] >= ORDINARY) || (key[pos1] >= ORDINARY))
	      {
		if(key[pos] == key[pos1])
		  {
		    key[pos] = EMPTY;
		    key[pos1] = EMPTY;
		    weight *= n_loop;
		    insert(1-data_position);
		  }
		else
		  {
		    join(pos,pos1);
		    key[pos] = EMPTY;
		    key[pos1] = EMPTY;
		    insert(1-data_position);
		  }
	      }
	    else if(key[pos] != key[pos1])
	      {
		key[pos] = EMPTY;
		key[pos1] = EMPTY;
		insert(1-data_position);
	      }

          temp=p_entry->p_next;
          free(p_entry);
          p_entry = temp;
        }
      p_hash[data_position][hash_entry]=NULL;
    }
  data_position = 1-data_position;
  
#ifdef DEBUG
  printf("after remove_aux_space : ");
  write_hash(data_position);
#endif
}

void remove_aux_space_special()
{
  int hash_entry,i;
  int pos,pos1;
  double loopfactor;
  T_entry* p_entry;
  T_entry* temp;

#ifdef DEBUG
  printf("before remove_aux_space_special : ");
  write_hash(data_position);
#endif

  pos = n;
  pos1 = n+1;
  
  for(hash_entry=0; hash_entry<hash_size; ++hash_entry)
    {
      p_entry = p_hash[data_position][hash_entry];
      while(p_entry)
        {
          unpack_state(p_entry);

          if((key[pos] == EMPTY) && (key[pos1] == EMPTY))
            insert_special(1-data_position);
          else if((key[pos] != EMPTY) && (key[pos1] != EMPTY))
	    if((key[pos]%32 >= ORDINARY) || (key[pos1]%32 >= ORDINARY))
	      {
		if(key[pos]%32 == key[pos1]%32)
		  {
		    loopfactor = getloopweight(key[pos]);
		    
		    key[pos] = EMPTY;
		    key[pos1] = EMPTY;
		    weight *= loopfactor;
		    insert_special(1-data_position);
		  }
		else
		  {
		    join(pos,pos1);
		    key[pos] = EMPTY;
		    key[pos1] = EMPTY;
		    insert_special(1-data_position);
		  }
	      }
	  /* Connect strings if they have matching target codes and stay clear of seam */
	    else if(((target[key[pos]%32]==key[pos1]%32) && (target[key[pos1]%32]==key[pos]%32)) &&
		    ((key[pos] >> 5) == (key[pos1] >> 5)))
	      {
		key[pos] = EMPTY;
		key[pos1] = EMPTY;
		insert_special(1-data_position);
	      }
	  
          temp=p_entry->p_next;
          free(p_entry);
          p_entry = temp;
        }
      p_hash[data_position][hash_entry]=NULL;
    }
  data_position = 1-data_position;
  
#ifdef DEBUG
  printf("after remove_aux_space_special : ");
  write_hash(data_position);
#endif
}

void transfer_row()
{
  int i;

  insert_aux_space();
  for(i=0; i<n; ++i)
    R_matrix(i);
  remove_aux_space();
  normalise_vector(data_position);
}

void transfer_row_special(int seam)
{
  int i;

  insert_aux_space_special(seam);
  for(i=0; i<n; ++i)
    R_matrix_special(i);
  remove_aux_space_special();
  normalise_vector(data_position);
}

void uncolor_strings()
{
  int hash_entry,i;
  int pos,pos1;
  T_entry* p_entry;
  T_entry* temp;

#ifdef DEBUG
  printf("before uncolor_strings : ");
  write_hash(data_position);
#endif

  for(hash_entry=0; hash_entry<hash_size; ++hash_entry)
    {
      p_entry = p_hash[data_position][hash_entry];
      while(p_entry)
        {
          unpack_state(p_entry);

	  for(i=0; i<n; ++i)
	    if(key[i] == MIDDLE)
	      key[i] = BOTTOM;
	  insert(1-data_position);

          temp=p_entry->p_next;
          free(p_entry);
          p_entry = temp;
        }
      p_hash[data_position][hash_entry]=NULL;
    }
  data_position = 1-data_position;
  
#ifdef DEBUG
  printf("after uncolor_strings : ");
  write_hash(data_position);
#endif
}

void uncolor_strings_special()
{
  int hash_entry,i;
  int pos,pos1;
  T_entry* p_entry;
  T_entry* temp;

#ifdef DEBUG
  printf("before uncolor_strings_special : ");
  write_hash(data_position);
#endif

  for(hash_entry=0; hash_entry<hash_size; ++hash_entry)
    {
      p_entry = p_hash[data_position][hash_entry];
      while(p_entry)
        {
          unpack_state(p_entry);

	  for(i=0; i<n; ++i)
	    if((key[i]%32 > EMPTY) && (key[i]%32 < ORDINARY))
	      key[i] = BOTTOM;
	  insert_special(1-data_position);

          temp=p_entry->p_next;
          free(p_entry);
          p_entry = temp;
        }
      p_hash[data_position][hash_entry]=NULL;
    }
  data_position = 1-data_position;
  
#ifdef DEBUG
  printf("after uncolor_strings_special : ");
  write_hash(data_position);
#endif
}

double compute_Z(int k1,int k2,int k3)
{
  int i;
  double result;

  /* Make sure that any inserted state has at least the required
     number of strings for the fusion to take place. In the bottom
     half this is k1 (indeed there is precisely k1 strings that
     cannot contract among themselves). In the top half, we can
     have between |k1-k2| and k1+k2 strings (in steps of two),
     but we must make sure that at least k3 strings survive to
     the end. If we don't do so, and if |k1-k2| < k3, then the
     states with k3 strings will decay exponentially and disappear
     from the state vector (due to finite numerical precision)! */

  /* Bottom operator */
  data_position = 0;
  min_strings = k1;
  for(i=0; i<k1; ++i)
    key[i] = BOTTOM;
  for(i=k1; i<N; ++i)
    key[i] = EMPTY;
  weight = 1.0;
  insert(data_position);
#ifdef VERBOSE
  write_hash(data_position);
#endif

  /* Bottom half cylinder */
  exponent = 0.0;
  for(i=0; i<M; ++i)
    {
      transfer_row();
#ifdef VERBOSE
      write_hash_statistics(data_position);
#endif
    }

  /* Middle operator */
  min_strings = k3;

  /* Top half cylinder */
 for(i=0; i<M; ++i)
    {
      transfer_row();
#ifdef VERBOSE
      write_hash_statistics(data_position);
#endif
    }

 /* Top operator */
 uncolor_strings();
 for(i=0; i<k3; ++i)
   key[i] = BOTTOM;
 for(i=k3; i<N; ++i)
   key[i] = EMPTY;
 result = get_weight(data_position);
 empty_hash(data_position);
 return result;
}

double compute_Z_special(int k1,int k2,int k3)
{
  int i,pairs;
  double result;

  /*
  printf("n=%.16f  n2=%.16f \n",n_loop,n2);
  printf("%d %d %d \n",k1,k2,k3);
  for(i=1;i<=9;++i)
    printf("%d %.16f \n",i,rho[i]);
  */
  
  /* Target codes */
  target[0] = -1; /* Not used */
  pairs = (k1-k3)/2;
  for(i=1; i<=k3; ++i)
    target[pairs+i] = k1+1; /* Bottom to top */
  for(i=0; i<pairs; ++i) /* Looping strings */
    {
      target[1+i] = k1-i;
      target[k1-i] = 1+i;
    }
  
  /* Bottom operator */
  data_position = 0;
  for(i=0; i<k1; ++i)
    key[i] = 1+i;
  for(i=k1; i<N; ++i)
    key[i] = EMPTY;
  weight = 1.0;
  insert_special(data_position);
#ifdef VERBOSE
  write_hash(data_position);
#endif
  
  /* Bottom half cylinder */
  exponent = 0.0;
  for(i=0; i<M; ++i)
    {
      transfer_row_special(1); /* Seam on */
#ifdef VERBOSE
      write_hash_statistics(data_position);
#endif
      // printf(" %.16f \n",exponent);
    }

  /* Middle operator */

  /* Top half cylinder */
 for(i=0; i<M; ++i)
    {
      transfer_row_special(0); /* Seam off */
#ifdef VERBOSE
      write_hash_statistics(data_position);
#endif
      // printf(" %.16f \n",exponent);
    }

 /* Top operator */
 uncolor_strings_special();
 for(i=0; i<k3; ++i)
   key[i] = BOTTOM;
 for(i=k3; i<N; ++i)
   key[i] = EMPTY;
 result = get_weight(data_position);
 empty_hash(data_position);
 return result;
}

double compute_Z_seam(int lower_seam,int upper_seam)
{
  int i;
  double result;

  /* Bottom operator */
  data_position = 0;
  for(i=0; i<N; ++i)
    key[i] = EMPTY;
  weight = 1.0;
  insert_special(data_position);
#ifdef VERBOSE
  write_hash(data_position);
#endif

  /* Bottom half cylinder */
  exponent = 0.0;
  for(i=0; i<M; ++i)
    {
      transfer_row_special(lower_seam);
#ifdef VERBOSE
      write_hash_statistics(data_position);
#endif
    }

  /* Middle operator */

  /* Top half cylinder */
 for(i=0; i<M; ++i)
    {
      transfer_row_special(upper_seam); /* Seam off */
#ifdef VERBOSE
      write_hash_statistics(data_position);
#endif
    }

 /* Top operator */
 uncolor_strings_special();
 for(i=0; i<N; ++i)
   key[i] = EMPTY;
 result = get_weight(data_position);
 empty_hash(data_position);
 return result;
}

#endif // TRANSFER_H_
