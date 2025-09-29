#include "common.h"
#include <stdlib.h>

double vector_norm()
{
  int i;
  double square_sum;
  T_entry* p_entry;
  
  square_sum = 0.0;
  for(i=0; i<hash_size; ++i)
    {
      p_entry = p_hash[data_position][i];
      while(p_entry)
	{
	  square_sum += (p_entry->weight) * (p_entry->weight);
	  p_entry = p_entry->p_next;
	}
    }
  return sqrt(square_sum);
}

void normalise_vector()
{
  int i;
  double v_norm;
  T_entry* p_entry;

  v_norm = vector_norm();
  exponent += log(v_norm);
  for(i=0; i<hash_size; ++i)
    {
      p_entry = p_hash[data_position][i];
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

void join(int i1,int i2)
{
  int i,val1,val2;
  
  if(key[i1] < key[i2])
    {
      val1 = key[i1];
      val2 = key[i2];
    }
  else
    {
      val1 = key[i2];
      val2 = key[i1];
    }

  for(i=0; i<N; ++i)
    if(key[i] == val2)
      key[i] = val1;
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

void insert_aux_space_special()
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
  double save_weight;
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
	  
          if((key[pos] == EMPTY) && (key[pos1] == EMPTY))
            {
              /* rho(1) delta(a,b,c,d) */
	      weight = save_weight * rho[1];
              insert_special(1-data_position); /* CORRECTED */

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
	      if((key[pos] >= ORDINARY) || (key[pos1] >= ORDINARY))
		{
		  if(key[pos] == key[pos1])
		    {
		      /* rho(4) delta(b,c,d) adj(a,b) */
		      key[pos] = EMPTY;
		      key[pos1] = EMPTY;
		      weight = save_weight * rho[4] * n_loop;
		      insert_special(1-data_position);
                  
		      /* rho(9) delta(b,d) adj(b,c) adj(a,b) */
		      key[pos] = ORDINARY+n+1;
		      key[pos1] = ORDINARY+n+1;
		      weight = save_weight * rho[9] * n_loop;
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
	      else if(((target[key[pos]]==0) && (target[key[pos1]]==kk1+2)) ||
		      ((target[key[pos]]==kk1+2) && (target[key[pos1]]==0)) ||
		      ((target[key[pos]]==key[pos1]) && (target[key[pos1]]==key[pos])))
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
	    if((key[pos] >= ORDINARY) || (key[pos1] >= ORDINARY))
	      {
		if(key[pos] == key[pos1])
		  {
		    key[pos] = EMPTY;
		    key[pos1] = EMPTY;
		    weight *= n_loop;
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
	    else if(((target[key[pos]]==0) && (target[key[pos1]]==kk1+2)) ||
		    ((target[key[pos]]==kk1+2) && (target[key[pos1]]==0)) ||
		    ((target[key[pos]]==key[pos1]) && (target[key[pos1]]==key[pos])))
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
  normalise_vector();
}

void transfer_row_special()
{
  int i;

  insert_aux_space_special();
  for(i=0; i<n; ++i)
    R_matrix_special(i);
  remove_aux_space_special();
  normalise_vector();
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
	    if((key[i] > EMPTY) && (key[i] < ORDINARY))
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
  printf("after uncolor_strings_special : ");
  write_hash(data_position);
#endif
}

int partner(int pos)
{
  int i;

  for(i=0; i<N; ++i)
    if((i != pos) && (key[i] == key[pos]))
      return i;
  printf("Cannot find partner to %d.\n",pos);
  exit(1);
}

void break_lines(int no_points)
{
  int hash_entry,i;
  int allowed,allowed0,pi;
  T_entry* p_entry;
  T_entry* temp;

#ifdef DEBUG
  printf("before break_lines(%d) : ",no_points);
  write_hash(data_position);
#endif

  for(hash_entry=0; hash_entry<hash_size; ++hash_entry)
    {
      p_entry = p_hash[data_position][hash_entry];
      while(p_entry)
        {
          unpack_state(p_entry);

	  // printf("treating "); write(key);

	  /* Break int(no_points/2) lines, by turning them into pairs of
	     dangling ends. Moreover, if no_points is odd, insert an edge
	     which is empty on the bottom and occupied on the top, AND
	     the converse, since we need the operator to be Hermitean! */

	  if((no_points % 2) == 0) // EVEN CASE
	    {
	      allowed = 1;
	      for(i=0; i<no_points/2; ++i)
		{
		  if(key[i] == EMPTY)
		    allowed = 0;
		  else if(key[i] >= ORDINARY)
		    {
		      pi = partner(i);
		      if(pi < no_points/2)
			allowed = 0;
		      else
			key[pi] = MIDDLE;
		    }
		}
	      if(allowed == 1)
		{
		  for(i=0; i<no_points/2; ++i)
		    key[i] = MIDDLE;
		  insert(1-data_position);
		}
	    }
	  else // ODD CASE
	    {
	      allowed0 = 1;
	      for(i=0; i<no_points/2 + 1; ++i)
		if(key[i] == EMPTY)
		  allowed0 = 0;
	      if(allowed0 == 1)
		{
		  // Case 1: Extra dangling end on the bottom
		  allowed = 1;
		  for(i=0; i<no_points/2 + 1; ++i)
		    {
		      if(key[i] >= ORDINARY)
			{
			  pi = partner(i);
			  if(pi < no_points/2 + 1)
			    allowed = 0;
			  else
			    key[pi] = MIDDLE;
			}
		    }
		  if(allowed == 1)
		    {
		      for(i=0; i<no_points/2; ++i)
			key[i] = MIDDLE;
		      key[no_points/2] = EMPTY;
		      insert(1-data_position);
		    }
		}
	      else
		{
		  allowed0 = 1;
		  for(i=0; i<no_points/2; ++i)
		    if(key[i] == EMPTY)
		      allowed0 = 0;
		  if(key[no_points/2] != EMPTY)
		    allowed0 = 0;
		  if(allowed0 == 1)
		    {
		      // Case 2: Extra dangling end on the top
		      allowed = 1;
		      for(i=0; i<no_points/2; ++i)
			{
			  if(key[i] >= ORDINARY)
			    {
			      pi = partner(i);
			      if(pi < no_points/2)
				allowed = 0;
			      else
				key[pi] = MIDDLE;
			    }
			}
		      if(allowed == 1)
			{
			  for(i=0; i<no_points/2 + 1; ++i)
			    key[i] = MIDDLE;
			  insert(1-data_position);
			}
		    }
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
  printf("after break_lines(%d) : ",no_points);  
  write_hash(data_position);
#endif
}

void break_lines_special(int no_points)
{
  int hash_entry,i;
  int allowed,allowed0,pi;
  T_entry* p_entry;
  T_entry* temp;

#ifdef DEBUG
  printf("before break_lines_special(%d) : ",no_points);
  write_hash(data_position);
#endif

  for(hash_entry=0; hash_entry<hash_size; ++hash_entry)
    {
      p_entry = p_hash[data_position][hash_entry];
      while(p_entry)
        {
          unpack_state(p_entry);

	  // printf("treating "); write(key);

	  /* Break int(no_points/2) lines, by turning them into pairs of
	     dangling ends. Moreover, if no_points is odd, insert an edge
	     which is empty on the bottom and occupied on the top, AND
	     the converse, since we need the operator to be Hermitean! */

	  if((no_points % 2) == 0) // EVEN CASE
	    {
	      allowed = 1;
	      for(i=0; i<no_points/2; ++i)
		{
		  if(key[i] == EMPTY)
		    allowed = 0;
		  else if(key[i] >= ORDINARY)
		    {
		      pi = partner(i);
		      if(pi < no_points/2)
			allowed = 0;
		      else
			key[pi] = kk1+1+i;
		    }
		  else if(target[key[i]] != 0)
		    allowed = 0;
		}
	      if(allowed == 1)
		{
		  for(i=0; i<no_points/2; ++i)
		    key[i] = kk1+no_points/2+1+i;
		  insert_special(1-data_position);
		}
	    }
	  else // ODD CASE
	    {
	      allowed0 = 1;
	      for(i=0; i<no_points/2 + 1; ++i)
		if(key[i] == EMPTY)
		  allowed0 = 0;
	      if(allowed0 == 1)
		{
		  // Case 1: Extra dangling end on the bottom
		  allowed = 1;
		  for(i=0; i<no_points/2 + 1; ++i)
		    {
		      if(key[i] >= ORDINARY)
			{
			  pi = partner(i);
			  if(pi < no_points/2 + 1)
			    allowed = 0;
			  else
			    key[pi] = kk1+1+i;
			}
		      else if(target[key[i]] != 0)
			allowed = 0;
		    }
		  if(allowed == 1)
		    {
		      for(i=0; i<no_points/2; ++i)
			key[i] = kk1+no_points/2+2+i;
		      key[no_points/2] = EMPTY;
		      insert_special(1-data_position);
		    }
		}
	      else
		{
		  allowed0 = 1;
		  for(i=0; i<no_points/2; ++i)
		    if(key[i] == EMPTY)
		      allowed0 = 0;
		  if(key[no_points/2] != EMPTY)
		    allowed0 = 0;
		  if(allowed0 == 1)
		    {
		      // Case 2: Extra dangling end on the top
		      allowed = 1;
		      for(i=0; i<no_points/2; ++i)
			{
			  if(key[i] >= ORDINARY)
			    {
			      pi = partner(i);
			      if(pi < no_points/2)
				allowed = 0;
			      else
				key[pi] = kk1+1+i;
			    }
			  else if(target[key[i]] != 0)
			    allowed = 0;
			}
		      if(allowed == 1)
			{
			  for(i=0; i<no_points/2 + 1; ++i)
			    key[i] = kk1+no_points/2+1+i;
			  insert_special(1-data_position);
			}
		    }
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
  printf("after break_lines_special(%d) : ",no_points);  
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
  min_strings = k1;
  for(i=0; i<k1; ++i)
    key[i] = BOTTOM;
  for(i=k1; i<n; ++i)
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
  break_lines(k2);
#ifdef VERBOSE
  write_hash_statistics(data_position);
#endif

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
 for(i=k3; i<n; ++i)
   key[i] = EMPTY;
 result = get_weight(data_position);
 empty_hash(data_position);
 return result;
}

double compute_Z_special(int k1,int k2,int k3)
{
  int i;
  double result;

  /* Target codes */
  target[0] = -1; /* Not used */
  for(i=1; i<=k1+k2; ++i)
    target[i] = 0; /* Bottom to middle */
  for(i=1; i<=k3; ++i)
    target[i] = k1+1; /* Bottom to top */
  for(i=0; i<(k1-k2-k3)/2; ++i) /* Looping strings */
    {
      target[k3+1+i] = k1-i;
      target[k1-i] = k3+1+i;
    }
  for(i=1; i<=k2; ++i)
    target[k1+i] = k1+2; /* Middle to bottom */

  /* Bottom operator */
  for(i=0; i<k1; ++i)
    key[i] = 1+i;
  for(i=k1; i<n; ++i)
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
      transfer_row_special();
#ifdef VERBOSE
      write_hash_statistics(data_position);
#endif
    }

  /* Middle operator */
  break_lines_special(k2);
#ifdef VERBOSE
  write_hash_statistics(data_position);
#endif

  /* Top half cylinder */
 for(i=0; i<M; ++i)
    {
      transfer_row_special();
#ifdef VERBOSE
      write_hash_statistics(data_position);
#endif
    }

 /* Top operator */
 uncolor_strings_special();
 for(i=0; i<k3; ++i)
   key[i] = BOTTOM;
 for(i=k3; i<n; ++i)
   key[i] = EMPTY;
 result = get_weight(data_position);
 empty_hash(data_position);
 return result;
}
