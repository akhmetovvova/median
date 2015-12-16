/*
 * ccarray.h
 *
 *  Created on: Dec 15, 2011
 *      Author: vova
 */

#ifndef __ccarray_h__
#define __ccarray_h__

#include <inttypes.h>
#include <stdlib.h>
#include <string.h>

#ifndef __symbian__
# include <malloc.h>
#endif



#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
  size_t capacity;
  size_t item_size;
  size_t items_count;
  void * items;
} ccarray_t;

typedef int (*cmpfunc_t)( const void * v1, const void * v2 );

static inline ccarray_t * ccarray_create( size_t capacity, size_t item_size )
{
  ccarray_t * c = (ccarray_t * )calloc(1, sizeof(ccarray_t));
  if ( c )
  {
    if ( !( c->items = calloc(capacity, item_size) ) )
    {
      free(c), c = 0;
    }
    else
    {
      c->capacity = capacity;
      c->items_count = 0;
      c->item_size = item_size;
    }
  }

  return c;
}

static inline void ccarray_destroy( ccarray_t * c )
{
  if ( c )
  {
    if ( c->items )
    {
      free(c->items);
    }

    free(c);
  }
}

static inline size_t ccarray_size( const ccarray_t * c )
{
  return c->items_count;
}

static inline size_t ccarray_capacity( const ccarray_t * c )
{
  return c->capacity;
}

static inline size_t ccarray_push_back( ccarray_t * c, const void * data )
{
  if ( c->items_count < c->capacity )
  {
    size_t pos = c->items_count++;
    memcpy((uint8_t*)c->items + pos * c->item_size, data, c->item_size);
    return pos;
  }
  return (size_t) ( -1 );
}

static inline size_t ccarray_push_front( ccarray_t * c, const void * data )
{
  if ( c->items_count < c->capacity )
  {
    memmove((uint8_t*)c->items + c->item_size, c->items, c->items_count * c->item_size);
    memcpy(c->items, data, c->item_size);
    return 0;
  }
  return (size_t) ( -1 );
}

static inline size_t ccarray_insert( ccarray_t * c, size_t pos, const void * data )
{
  if ( c->items_count < c->capacity )
  {
    memmove((uint8_t*)c->items + ( pos + 1 ) * c->item_size, (uint8_t*)c->items + pos * c->item_size,
      ( c->items_count++ - pos ) * c->item_size);
    memcpy((uint8_t*)c->items + pos * c->item_size, data, c->item_size);
    return pos;
  }
  return (size_t) ( -1 );
}

static inline size_t ccarray_pop_back( ccarray_t * c, void * data )
{
  if ( c->items_count < 1 )
  {
    return (size_t) ( -1 );
  }

  memcpy(data, (uint8_t*)c->items + --c->items_count * c->item_size, c->item_size);

  return c->items_count;
}

static inline size_t ccarray_pop_front( ccarray_t * c, void * data )
{
  if ( c->items_count < 1 ) {
    return (size_t) ( -1 );
  }

  if ( data ) {
    memcpy(data, c->items, c->item_size);
  }

  memcpy(c->items, (uint8_t*)c->items + c->item_size, --c->items_count * c->item_size);
  return c->items_count;
}

static inline size_t ccarray_erase( ccarray_t * c, size_t pos )
{
  if ( pos < c->items_count )
  {
    memcpy((uint8_t*)c->items + pos * c->item_size, (uint8_t*)c->items + ( pos + 1 ) * c->item_size,
      ( --c->items_count - pos ) * c->item_size);
  }

  return c->items_count;
}

static inline void ccarray_clear( ccarray_t * c )
{
  c->items_count = 0;
}

static inline void * ccarray_peek( const ccarray_t * c, size_t pos )
{
  return ((uint8_t*)c->items) + pos * c->item_size;
}

static inline void * ccarray_ppeek( ccarray_t * c, size_t pos )
{
  return *(void **) ( (uint8_t*) c->items + pos * c->item_size );
}

static inline void * ccarray_peek_end( ccarray_t * c )
{
  return (uint8_t*)c->items + c->items_count * c->item_size;
}

static inline size_t ccarray_set_size( ccarray_t * c, size_t newsize )
{
  return (c->items_count = newsize);
}


static inline size_t ccarray_find_item( ccarray_t * c, const void * value )
{
  size_t pos = 0;
  const uint8_t * begin = (uint8_t*)c->items;
  const uint8_t * const end = (uint8_t*)c->items + c->items_count * c->item_size;
  while ( begin < end && bcmp(begin, value, c->item_size) != 0 )
  {
    begin += c->item_size, ++pos;
  }
  return pos;
}

static inline size_t ccarray_find( ccarray_t * c, cmpfunc_t cmp, const void * value )
{
  size_t pos = 0;
  const uint8_t * begin = (uint8_t*)c->items;
  const uint8_t * const end = (uint8_t*)c->items + c->items_count * c->item_size;
  while ( begin < end && cmp(begin, value) != 0 )
  {
    begin += c->item_size, ++pos;
  }
  return pos;
}

static inline void ccarray_sort( ccarray_t * c, size_t beg, size_t end, cmpfunc_t cmp)
{
  qsort(((uint8_t*) c->items) + beg * c->item_size, end - beg, c->item_size, cmp);
}

static inline size_t ccarray_lowerbound(const ccarray_t * c, size_t beg, size_t end, cmpfunc_t cmp, const void * value)
{
  const size_t origbeg = beg;
  size_t len = end - beg;
  size_t half, mid;
  int rc;

  while ( len > 0 )
  {
    mid = beg + ( half = len >> 1 );

    rc = cmp((uint8_t*)c->items + mid * c->item_size, value);

    if ( rc == 0 ) {
      beg = mid;
      break;
    }

    if ( rc > 0 ) {
      len = half;
    }
    else {
      beg = mid + 1;
      len -= ( half + 1 );
    }
  }

  while ( beg > origbeg && cmp((uint8_t*) c->items + (beg - 1) * c->item_size, value) == 0 ) {
    -- beg;
  }

  return beg;
}

static inline int ccarray_create_slots( ccarray_t ** a, ccarray_t ** apool, size_t num_slots, size_t item_size )
{
  if ( !( *a = ccarray_create(num_slots, item_size) ) )
  {
    return -1;
  }

  if ( !( *apool = ccarray_create(num_slots, item_size) ) )
  {
    ccarray_destroy(*a), *a = 0;
    return -1;
  }

  ( *apool )->items_count = num_slots;

  return 0;
}

static inline void ccarray_destroy_slots( ccarray_t ** a, ccarray_t ** apool )
{
  if ( *a )
  {
    ccarray_destroy(*a);
    *a = 0;
  }

  if ( *apool )
  {
    ccarray_destroy(*apool);
    *apool = 0;
  }
}


static inline int ccarray_create_pslots( ccarray_t ** aq, ccarray_t ** ap, size_t num_slots, size_t slot_size )
{
  size_t i, j;

  if ( !( *aq = ccarray_create(num_slots, sizeof(void *)) ) )
  {
    return -1;
  }

  if ( !( *ap = ccarray_create(num_slots, sizeof(void *)) ) )
  {
    ccarray_destroy(*aq), *aq = 0;
    return -1;
  }

  for ( i = 0; i < num_slots; ++i )
  {
    void * p = calloc(1, slot_size);
    if ( p )
    {
      ccarray_push_back(*ap, &p);
    }
    else
    {
      for ( j = 0; j < i; ++j )
      {
        ccarray_pop_back(*ap, &p);
        free(p);
      }

      ccarray_destroy(*aq), *aq = 0;
      ccarray_destroy(*ap), *ap = 0;
      return -1;
    }
  }

  return 0;
}

static inline void ccarray_destroy_pslots( ccarray_t ** aq, ccarray_t ** ap )
{
  void * p = 0;

  if ( *aq )
  {
    while ( ccarray_size(*aq) > 0 )
    {
      ccarray_pop_back(*aq, &p);
      if ( p )
      {
        free(p);
      }
    }
    ccarray_destroy(*aq), *aq = 0;
  }

  if ( *ap )
  {
    while ( ccarray_size(*ap) > 0 )
    {
      ccarray_pop_back(*ap, &p);
      if ( p )
      {
        free(p);
      }
    }
    ccarray_destroy(*ap), *ap = 0;
  }
}


#ifdef __cplusplus
}
#endif

#endif /* __ccarray_h__ */
