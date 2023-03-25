#include "rope.h"

#include <iostream>
#include <vector>

#include "CGL/vector2D.h"
#include "mass.h"
#include "spring.h"

namespace CGL {

Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass,
           float k, vector<int> pinned_nodes) {
  // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and
  // containing `num_nodes` nodes.
  //
  // New masses can be created (see mass.h):
  // Mass *m = new Mass(position, mass, bool)
  // and should be added to the masses vector.
  //
  // Springs can be created (see spring.h):
  // Spring *s = new Spring(mass_a, mass_b, spring_const)
  // and should be added to the springs vector.
  //
  // Masses corresponding to indices in pinned_nodes
  // should have their pinned field set to true.

  Vector2D width = (start - end) / (num_nodes - 1);

  for (int i = 0; i <= num_nodes; i++) {
    Vector2D position = start + (i * width);
    Mass *mass = new Mass(position, node_mass, false);
    masses.push_back(mass);

    // Pair previous mass with current mass to create spring.
    if (i > 0) {
      Spring *spring = new Spring(masses[i - 1], mass, k);
      springs.push_back(spring);
    }
  }

  for (int i : pinned_nodes) {
    masses[i]->pinned = true;
  }
}

void Rope::simulateEuler(float delta_t, Vector2D gravity) {
  for (auto &s : springs) {
    // TODO (Part 2.1): Use Hooke's law to calculate the force on a node

    // s->m1 is a and s->m2 is b
    Vector2D dist = s->m2->position - s->m1->position;
    double dist_norm = dist.norm();
    // spring force f_a->b
    Vector2D f_a_b = s->k * (dist / dist_norm) * (dist_norm - s->rest_length);
    s->m1->forces += f_a_b;
    s->m2->forces -= f_a_b;

    // TODO (Part 4.1): Add damping forces

    // TODO Apply forces as appropriate.
  }

  for (auto &m : masses) {
    if (!m->pinned) {
      // TODO (Part 2.1): Add the force due to gravity, then compute the new
      // velocity and position

      // F = ma => a = F/m
      Vector2D acceleration = m->forces / m->mass + gravity;
      m->velocity = m->velocity + acceleration * delta_t;
      m->position = m->position + m->velocity * delta_t;
    }

    // TODO Reset all forces on each mass
    m->forces = Vector2D(0.0, 0.0);
  }
}

void Rope::simulateVerlet(float delta_t, Vector2D gravity) {
  // TODO (Part 3.1): Clear forces

  for (auto &m : masses) {
    m->forces = Vector2D(0.0, 0.0);
  }

  for (auto &s : springs) {
    // TODO (Part 3.1): Simulate one timestep of the rope using explicit Verlet

    // s->m1 is a and s->m2 is b
    Vector2D dist = s->m2->position - s->m1->position;
    double dist_norm = dist.norm();
    // spring force f_a->b
    Vector2D f_a_b = s->k * (dist / dist_norm) * (dist_norm - s->rest_length);
    s->m1->forces += f_a_b;
    s->m2->forces -= f_a_b;
  }

  for (auto &m : masses) {
    if (!m->pinned) {
      Vector2D temp_position = m->position;

      // TODO (Part 3.1): Set the new position of the rope mass

      // F = ma => a = F/m
      Vector2D acceleration = m->forces / m->mass + gravity;
      // m->position = m->position + (m->position - m->last_position) +
      //               acceleration * delta_t * delta_t;

      // TODO (Part 4.2): Add global Verlet damping

      double damping_factor = 0.0003;
      m->position = m->position +
                    (1 - damping_factor) * (m->position - m->last_position) +
                    acceleration * delta_t * delta_t;

      m->last_position = temp_position;
    }
  }
}
}  // namespace CGL
