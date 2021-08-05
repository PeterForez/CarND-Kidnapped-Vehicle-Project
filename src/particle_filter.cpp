/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>     // Need this for sampling from distributions
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) 
{
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   * @param x       GPS provided x position
   * @param y       GPS provided y position
   * @param theta   GPS provided yaw
   * @param std[]   Standard deviation of x, y, theta
   */
  num_particles = 100;                                     // Set the number of particles
  Particle particle;
  std::cout << "num_particles " << num_particles << std::endl;
  
  double std_x     = std[0];                               // Standard deviations for x
  double std_y     = std[1];                               // Standard deviations for y
  double std_theta = std[2];                               // Standard deviations for theta

  normal_distribution<double> dist_x(x, std_x);            // Create normal distributions for x
  normal_distribution<double> dist_y(y, std_y);            // Create normal distributions for y
  normal_distribution<double> dist_theta(theta, std_theta);// Create normal distributions for theta
  
  std::default_random_engine gen;                          // This is a random number engine class that generates pseudo-random numbers.
  
  for (int i = 0; i < num_particles; ++i) 
  {
    particle.id     = i;                                   // Id of the particle in the map.
    particle.x      = dist_x(gen);                         // Sample from the normal distribution of x
    particle.y      = dist_y(gen);                         // Sample from the normal distribution of y
    particle.theta  = dist_theta(gen);                     // Sample from the normal distribution of theta
    particle.weight = 1;                                   // Initialize all weights to 1
    
    particles.push_back(particle);                         // Vector of all the particles created
    
    //std::cout << "Sample " << i + 1 << " " << particle.x << " " << particle.y << " " << particle.theta << std::endl;
  }
  
  is_initialized = true;
  return;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) 
{
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  double std_x     = std[0];                               // Standard deviations for x
  double std_y     = std[1];                               // Standard deviations for y
  double std_theta = std[2];                               // Standard deviations for theta
  
  normal_distribution<double> dist_x(x, std_x);            // Create normal distributions for x
  normal_distribution<double> dist_y(y, std_y);            // Create normal distributions for y
  normal_distribution<double> dist_theta(theta, std_theta);// Create normal distributions for theta
  
  std::default_random_engine gen;                          // This is a random number engine class that generates pseudo-random numbers.
    
  for (int i; i < particles.size(); i++)
  {
    if(fabs(yaw_rate) > 0.0001)                            // Absolute yaw rate is not equal to zero
    {
      particles[i].x     += velocity / yaw_rate * [sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta)];
      particles[i].y     += velocity / yaw_rate * [cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t)];
      particles[i].theta += yaw_rate * delta_t;
    }
    else                                                   // Yaw rate is equal to zero
    {
      particles[i].x     += velocity * delta_t * cos(particles[i].theta);
      particles[i].y     += velocity * delta_t * sin(particles[i].theta);
    }
    // Adding Noise
    particles[i].x       += dist_x(gen);                   // Sample from the normal distribution of x
    particles[i].y       += dist_y(gen);                   // Sample from the normal distribution of y
    particles[i].theta   += dist_theta(gen);               // Sample from the normal distribution of theta
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) 
{
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
   
  //https://knowledge.udacity.com/questions/516274 
  for (unsigned int i; i < observations.size(); i++)                                         // Loop over the Observations
  {                                                                                          
    double min_distance = numeric_limits<double>::max();                                     // Initialize with the maximum value of double
    for (unsigned int j; j < predicted.size(); j++)                                          // Loop over the Predictions    
    {                                                                                        
      double distance;                                                                       
      distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y); // Function in "helper_functions.h"
      observations[i].id = -1;                                                               // Initialize the observation id
      if (distance < min_distance)                                                           // Check the minimum distance
      {                                                                                      
        min_distance = distance;                                                             
        observations[i].id = predicted[j].id;                                                // Associating the observation to landmark id
      }
    }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) 
{
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  
}

void ParticleFilter::resample() 
{
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) 
{
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) 
{
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) 
{
  vector<double> v;

  if (coord == "X") 
  {
    v = best.sense_x;
  } 
  else 
  {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}