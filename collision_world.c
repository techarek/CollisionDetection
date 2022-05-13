/**
 * collision_world.c -- detect and handle line segment intersections
 * Copyright (c) 2012-2020 the Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 **/

#include "./collision_world.h"

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <cilk/cilk.h>
#include <cilk/cilkscale.h>
#include <cilk/reducer.h>

#include "./intersection_detection.h"
#include "./intersection_event_list.h"
#include "./line.h"

typedef CILK_C_DECLARE_REDUCER(IntersectionEventList) IntersectionEventListReducer;

#define NUM_STRIPS 8
#define EPS 1e-6
#define COARSEN 1000
#define MAX(a, b) (((a) > (b))? (a): (b))
#define MIN(a, b) (((a) < (b))? (a): (b))


void my_qsort(Line** begin, Line** end) {
  if (begin >= end - 1) {
    return;
  }

  // Perform a partition on the first element of the array
  // so that the left consists of elements less than 'begin'
  // and the right consists of elements greater than 'begin'
  Line** scan = begin + 1;
  Line** swap = end;

  while (scan < swap) {
    if ((*scan)->min_x > (*begin)->min_x) {
      swap--;
      Line* tmp = *swap;
      *swap = *scan;
      *scan = tmp;
    } else {
      scan++;
    }
  }

  Line** middle = swap - 1;  // Last element of left partition.

  // Swap the element we partitioned on into the middle.
  Line* tmp = *begin;
  *begin = *middle;
  *middle = tmp;

  Line** maximum = MAX(begin + 1, middle);

  // Recurse on the left and right partitions.
  cilk_spawn my_qsort(begin, middle);
  my_qsort(maximum, end);
  cilk_sync;
}

CollisionWorld* CollisionWorld_new(const unsigned int capacity) {
  assert(capacity > 0);

  CollisionWorld* collisionWorld = malloc(sizeof(CollisionWorld));
  if (collisionWorld == NULL) {
    return NULL;
  }

  collisionWorld->numLineWallCollisions = 0;
  collisionWorld->numLineLineCollisions = 0;
  collisionWorld->timeStep = 0.5;
  collisionWorld->lines = malloc(capacity * sizeof(Line*));
  collisionWorld->numOfLines = 0;
  return collisionWorld;
}

void CollisionWorld_delete(CollisionWorld* collisionWorld) {
  #pragma cilk grainsize COARSEN
  cilk_for (int i = 0; i < collisionWorld->numOfLines; i++) {
    free(collisionWorld->lines[i]);
  }
  cilk_sync;
  free(collisionWorld->lines);
  free(collisionWorld);
}

unsigned int CollisionWorld_getNumOfLines(CollisionWorld* collisionWorld) {
  return collisionWorld->numOfLines;
}

void CollisionWorld_addLine(CollisionWorld* collisionWorld, Line *line) {
  collisionWorld->lines[collisionWorld->numOfLines] = line;
  collisionWorld->numOfLines++;
}

Line* CollisionWorld_getLine(CollisionWorld* collisionWorld,
                             const unsigned int index) {
  if (index >= collisionWorld->numOfLines) {
    return NULL;
  }
  return collisionWorld->lines[index];
}

void CollisionWorld_updateLines(CollisionWorld* collisionWorld) {
  CollisionWorld_detectIntersection(collisionWorld);
  CollisionWorld_updatePosition(collisionWorld);
  CollisionWorld_lineWallCollision(collisionWorld);
}

void CollisionWorld_updatePosition(CollisionWorld* collisionWorld) {
  double t = collisionWorld->timeStep;
  #pragma cilk grainsize COARSEN
  cilk_for (int i = 0; i < collisionWorld->numOfLines; i++) {
    Line *line = collisionWorld->lines[i];
    line->p1 = Vec_add(line->p1, Vec_multiply(line->velocity, t));
    line->p2 = Vec_add(line->p2, Vec_multiply(line->velocity, t));
  }
  cilk_sync;
}

void CollisionWorld_lineWallCollision(CollisionWorld* collisionWorld) {
  for (int i = 0; i < collisionWorld->numOfLines; i++) {
    Line *line = collisionWorld->lines[i];
    bool collide = false;

    // Right side
    if ((line->p1.x > BOX_XMAX || line->p2.x > BOX_XMAX)
        && (line->velocity.x > 0)) {
      line->velocity.x = -line->velocity.x;
      collide = true;
    }
    // Left side
    if ((line->p1.x < BOX_XMIN || line->p2.x < BOX_XMIN)
        && (line->velocity.x < 0)) {
      line->velocity.x = -line->velocity.x;
      collide = true;
    }
    // Top side
    if ((line->p1.y > BOX_YMAX || line->p2.y > BOX_YMAX)
        && (line->velocity.y > 0)) {
      line->velocity.y = -line->velocity.y;
      collide = true;
    }
    // Bottom side
    if ((line->p1.y < BOX_YMIN || line->p2.y < BOX_YMIN)
        && (line->velocity.y < 0)) {
      line->velocity.y = -line->velocity.y;
      collide = true;
    }
    // Update total number of collisions.
    if (collide == true) {
       collisionWorld->numLineWallCollisions++;
    }
  }
}

void sweepLine(CollisionWorld* collisionWorld, IntersectionEventListReducer* reducer_ptr, int strip) {

  double strip_start = (strip == 0)? 0:BOX_YMIN + 0.5/NUM_STRIPS * strip;
  double strip_end = (strip == NUM_STRIPS - 1)? 2:BOX_YMIN + 0.5/NUM_STRIPS * (strip + 1);

  Line* in_strip[collisionWorld->numOfLines];
  int num_in_strip = 0;

  for (int i = 0; i < collisionWorld->numOfLines; i++) {
    Line* line = collisionWorld->lines[i];
    if (!(line->min_y > strip_end || strip_start > line->max_y)) {
      in_strip[num_in_strip++] = line;
    }
  }

  my_qsort(in_strip, in_strip + num_in_strip);

  // SWEEP LINE ALGORITHM
  cilk_for (int i = 0; i < num_in_strip; i++) {
    for (int j = i + 1; j < num_in_strip; j++) {
      Line *l1 = in_strip[i];
      Line *l2 = in_strip[j];

        if (l2->min_x > l1->max_x) {
          break;
        } 

        // intersect expects compareLines(l1, l2) < 0 to be true.
        // Swap l1 and l2, if necessary.
        if (compareLines(l1, l2) >= 0) {
          Line *temp = l1;
          l1 = l2;
          l2 = temp;
        }

        IntersectionType intersectionType =
            intersect(l1, l2, collisionWorld->timeStep);
        if (intersectionType != NO_INTERSECTION) {
          IntersectionEventList_appendNode(&REDUCER_VIEW(*reducer_ptr), l1, l2,
                                          intersectionType);
        }
    }
  }
  cilk_sync;
}

IntersectionEventListReducer X = CILK_C_INIT_REDUCER(IntersectionEventList,
IntersectionEventList_reduce, IntersectionEventList_identity, IntersectionEventList_destroy,
(IntersectionEventList) { .head = NULL, .tail = NULL, .size = 0 }); 

void CollisionWorld_detectIntersection(CollisionWorld* collisionWorld) {

  double t = collisionWorld->timeStep;
  #pragma cilk grainsize COARSEN
  cilk_for (int i = 0; i < collisionWorld->numOfLines; ++i) {
    Line* line = collisionWorld->lines[i];
    Vec p3 = Vec_add(line->p1, Vec_multiply(line->velocity, t));
    Vec p4 = Vec_add(line->p2, Vec_multiply(line->velocity, t));
    
    line->min_x = MIN(MIN(line->p1.x, line->p2.x), MIN(p3.x, p4.x)) - EPS;
    line->max_x = MAX(MAX(line->p1.x, line->p2.x), MAX(p3.x, p4.x)) + EPS;
    line->min_y = MIN(MIN(line->p1.y, line->p2.y), MIN(p3.y, p4.y)) - EPS;
    line->max_y = MAX(MAX(line->p1.y, line->p2.y), MAX(p3.y, p4.y)) + EPS;
  }
  cilk_sync;

  CILK_C_REGISTER_REDUCER(X);

  cilk_for (int i = 0; i < NUM_STRIPS; ++i) {
    sweepLine(collisionWorld, &X, i);
  }

  cilk_sync;
 
  #ifndef NDEBUG
    IntersectionEventList intersectionEventListOg = IntersectionEventList_make();
    int numCollisions = 0;
    // Original Implementation
    for (int i = 0; i < collisionWorld->numOfLines; i++) {
      for (int j = i + 1; j < collisionWorld->numOfLines; j++) {
        Line *l1 = collisionWorld->lines[i];
        Line *l2 = collisionWorld->lines[j];

        // intersect expects compareLines(l1, l2) < 0 to be true.
        // Swap l1 and l2, if necessary.
        if (compareLines(l1, l2) >= 0) {
          Line *temp = l1;
          l1 = l2;
          l2 = temp;
        }

        IntersectionType intersectionType =
            intersect(l1, l2, collisionWorld->timeStep);
        if (intersectionType != NO_INTERSECTION) {
          IntersectionEventList_appendNode(&intersectionEventListOg, l1, l2,
                                          intersectionType);
          numCollisions++;
        }
      }
    }

  #endif


  // Sort the intersection event list.
  IntersectionEventNode* startNode = X.value.head;
  while (startNode != NULL) {
    IntersectionEventNode* minNode = startNode;
    IntersectionEventNode* curNode = startNode->next;
    while (curNode != NULL) {
      if (IntersectionEventNode_compareData(curNode, minNode) < 0) {
        minNode = curNode;
      }
      curNode = curNode->next;
    }
    if (minNode != startNode) {
      IntersectionEventNode_swapData(minNode, startNode);
    }
    startNode = startNode->next;
  }

  // Call the collision solver for each intersection event.
  IntersectionEventNode* curNode = X.value.head;

  while (curNode != NULL) {
    CollisionWorld_collisionSolver(collisionWorld, curNode->l1, curNode->l2,
                                  curNode->intersectionType);

    collisionWorld->numLineLineCollisions++;

    while (curNode->next != NULL && curNode->next->l1->id == curNode->l1->id && curNode->next->l2->id == curNode->l2->id) {
      curNode = curNode->next;
      X.value.size--;
    }
    curNode = curNode->next;

  }

  #ifndef NDEBUG
    printf("Expected: %d, Returned: %d\n", numCollisions, X.value.size);
    assert (X.value.size == numCollisions);
  #endif
  
  IntersectionEventList_deleteNodes(&X.value);
  CILK_C_UNREGISTER_REDUCER(X);
}

unsigned int CollisionWorld_getNumLineWallCollisions(
    CollisionWorld* collisionWorld) {
  return collisionWorld->numLineWallCollisions;
}

unsigned int CollisionWorld_getNumLineLineCollisions(
    CollisionWorld* collisionWorld) {
  return collisionWorld->numLineLineCollisions;
}

void CollisionWorld_collisionSolver(CollisionWorld* collisionWorld,
                                    Line *l1, Line *l2,
                                    IntersectionType intersectionType) {
  assert(compareLines(l1, l2) < 0);
  assert(intersectionType == L1_WITH_L2
         || intersectionType == L2_WITH_L1
         || intersectionType == ALREADY_INTERSECTED);

  // Despite our efforts to determine whether lines will intersect ahead
  // of time (and to modify their velocities appropriately), our
  // simplified model can sometimes cause lines to intersect.  In such a
  // case, we compute velocities so that the two lines can get unstuck in
  // the fastest possible way, while still conserving momentum and kinetic
  // energy.
  if (intersectionType == ALREADY_INTERSECTED) {
    Vec p = getIntersectionPoint(l1->p1, l1->p2, l2->p1, l2->p2);

    if (Vec_length(Vec_subtract(l1->p1, p))
        < Vec_length(Vec_subtract(l1->p2, p))) {
      l1->velocity = Vec_multiply(Vec_normalize(Vec_subtract(l1->p2, p)),
                                  Vec_length(l1->velocity));
    } else {
      l1->velocity = Vec_multiply(Vec_normalize(Vec_subtract(l1->p1, p)),
                                  Vec_length(l1->velocity));
    }
    if (Vec_length(Vec_subtract(l2->p1, p))
        < Vec_length(Vec_subtract(l2->p2, p))) {
      l2->velocity = Vec_multiply(Vec_normalize(Vec_subtract(l2->p2, p)),
                                  Vec_length(l2->velocity));
    } else {
      l2->velocity = Vec_multiply(Vec_normalize(Vec_subtract(l2->p1, p)),
                                  Vec_length(l2->velocity));
    }
    return;
  }

  // Compute the collision face/normal vectors.
  Vec face;
  Vec normal;
  if (intersectionType == L1_WITH_L2) {
    Vec v = Vec_makeFromLine(*l2);
    face = Vec_normalize(v);
  } else {
    Vec v = Vec_makeFromLine(*l1);
    face = Vec_normalize(v);
  }
  normal = Vec_orthogonal(face);

  // Obtain each line's velocity components with respect to the collision
  // face/normal vectors.
  double v1Face = Vec_dotProduct(l1->velocity, face);
  double v2Face = Vec_dotProduct(l2->velocity, face);
  double v1Normal = Vec_dotProduct(l1->velocity, normal);
  double v2Normal = Vec_dotProduct(l2->velocity, normal);

  // Compute the mass of each line (we simply use its length).
  double m1 = Vec_length(Vec_subtract(l1->p1, l1->p2));
  double m2 = Vec_length(Vec_subtract(l2->p1, l2->p2));

  // Perform the collision calculation (computes the new velocities along
  // the direction normal to the collision face such that momentum and
  // kinetic energy are conserved).
  double newV1Normal = ((m1 - m2) / (m1 + m2)) * v1Normal
      + (2 * m2 / (m1 + m2)) * v2Normal;
  double newV2Normal = (2 * m1 / (m1 + m2)) * v1Normal
      + ((m2 - m1) / (m2 + m1)) * v2Normal;

  // Combine the resulting velocities.
  l1->velocity = Vec_add(Vec_multiply(normal, newV1Normal),
                         Vec_multiply(face, v1Face));
  l2->velocity = Vec_add(Vec_multiply(normal, newV2Normal),
                         Vec_multiply(face, v2Face));

  return;
}
