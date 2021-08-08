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
  
  //std::cout << "Step: ParticleFilter::init " << std::endl;
  
  num_particles = 100;                                     // Set the number of particles
  std::cout << "num_particles " << num_particles << std::endl;
  
  double std_x     = std[0];                               // Standard deviations for x
  double std_y     = std[1];                               // Standard deviations for y
  double std_theta = std[2];                               // Standard deviations for theta

  std::normal_distribution<double> dist_x(x, std_x);            // Create normal distributions for x
  std::normal_distribution<double> dist_y(y, std_y);            // Create normal distributions for y
  std::normal_distribution<double> dist_theta(theta, std_theta);// Create normal distributions for theta
  
  std::default_random_engine gen;                          // This is a random number engine class that generates pseudo-random numbers.
  Particle particle;
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
  
  is_initialized = true;                                   // Flag to indicate that the step of initialization is finished
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
  
  //std::cout << "Step: ParticleFilter::prediction " << std::endl;
  
  double std_x     = std_pos[0];                           // Standard deviations for x
  double std_y     = std_pos[1];                           // Standard deviations for y
  double std_theta = std_pos[2];                           // Standard deviations for theta
  
  std::normal_distribution<double> dist_x(0, std_x);            // Create normal distributions for x with zero mean
  std::normal_distribution<double> dist_y(0, std_y);            // Create normal distributions for y with zero mean
  std::normal_distribution<double> dist_theta(0, std_theta);    // Create normal distributions for theta with zero mean
  
  std::default_random_engine gen;                          // This is a random number engine class that generates pseudo-random numbers.
    
  for (size_t i = 0; i < particles.size(); i++)
  {
    if(fabs(yaw_rate) > 0.0001)                            // Absolute yaw rate is not equal to zero
    {
      particles[i].x     += (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
      particles[i].y     += (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
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

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, vector<LandmarkObs>& observations) 
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
  
  //std::cout << "Step: ParticleFilter::dataAssociation " << std::endl;
  
  for (size_t i = 0; i < observations.size(); i++)                                           // Loop over the Observations
  {                                                                                          
    double min_distance = std::numeric_limits<double>::max();                                     // Initialize with the maximum value of double
    for (size_t j = 0; j < predicted.size(); j++)                                            // Loop over the Predictions    
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

// https://classroom.udacity.com/nanodegrees/nd013/parts/b9040951-b43f-4dd3-8b16-76e7b52f4d9d/modules/85ece059-1351-4599-bb2c-0095d6534c8c/lessons/e3981fd5-8266-43be-a497-a862af9187d4/concepts/0a756b5c-458b-491f-b560-ac18b251f14d

/**
 * x and y are the observations in map coordinates from landmarks and μx and μy
 * are the coordinates of the nearest landmarks. 
 */
double multiv_prob(double sig_x, double sig_y, 
                   double x_obs, double y_obs,
                   double mu_x,  double mu_y) 
{
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
           + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
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
   
  //std::cout << "Step: ParticleFilter::updateWeights " << std::endl;
  
  double x_obs;      // The x coordinate for the landmark observation
  double y_obs;      // The y coordinate for the landmark observation
                     
  double x_part;     // The x coordinate for the particle
  double y_part;     // The y coordinate for the particle
  double theta;      // The theta orientation of the particle
                     
  double x_map;      // The x transform to map coordinate for the landmark observation
  double y_map;      // The y transform to map coordinate for the landmark observation
  
  float x_landmark;  // Landmark x-position in the map (global coordinates)
  float y_landmark;  // Landmark y-position in the map (global coordinates)
  
  double mu_x;
  double mu_y;
  
  int   id ;         // Landmark or Observation ID
  
  double sig_x = std_landmark[0];
  double sig_y = std_landmark[1];
  
  double distance;  

  for (size_t i = 0; i < particles.size(); i++)
  {
    x_part = particles[i].x;
    y_part = particles[i].y;
    theta  = particles[i].theta; 
    
    // Find Landmarks within the sensor range of the particle
    vector<LandmarkObs> landmarks_sensor_range;
    for (size_t j = 1; j < map_landmarks.landmark_list.size(); j++)
    {
      id         = map_landmarks.landmark_list[j].id_i; 
      x_landmark = map_landmarks.landmark_list[j].x_f; 
      y_landmark = map_landmarks.landmark_list[j].y_f; 
      
      distance   = dist(x_part, x_landmark, y_part, y_landmark);
      if (distance <= sensor_range) // Landmark within the sensor range
      {
        landmarks_sensor_range.push_back(LandmarkObs{id, x_landmark, y_landmark}); // Prediction
      }
    }
    
    // Transform the observation from Vehicile to MAP coordinates    
    vector<LandmarkObs> transformed_observations;
    for (size_t j = 0; j < observations.size(); j++)
    {
      id    = observations[i].id;
      x_obs = observations[i].x;
      y_obs = observations[i].y;
      
      x_map = x_part + (cos(theta) * x_obs) - (sin(theta) * y_obs);
      y_map = y_part + (sin(theta) * x_obs) + (cos(theta) * y_obs);
      transformed_observations.push_back(LandmarkObs{id, x_map, y_map}); 
    }
    
    // Find the predicted measurement that is closest to each observed measurement and assign the observed measurement to this particular landmark.
    dataAssociation(landmarks_sensor_range, transformed_observations);
    
    // Update the Weight with multivaraiate gaussian distribution
    particles[i].weight = 1.0; // Initialize the Weights to one
    for (size_t j = 0; j < transformed_observations.size(); j++)
    {
      x_obs = transformed_observations[j].x;
      y_obs = transformed_observations[j].y;
      id    = transformed_observations[j].id;
      for (size_t k = 0; k < landmarks_sensor_range.size(); k++) 
      {
        if (landmarks_sensor_range[k].id == id) 
        {
          mu_x = landmarks_sensor_range[k].x;
          mu_y = landmarks_sensor_range[k].y;
          break;
        }
      }
      particles[i].weight *= multiv_prob(sig_x, sig_y, x_obs, y_obs, mu_x, mu_y);
    }
    
    
  }
}

void ParticleFilter::resample() 
{
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
   
  // https://knowledge.udacity.com/questions/240067
  
  //std::cout << "Step: ParticleFilter::resample " << std::endl;
  
  // Maximum and Total Weight
  double weight_max = std::numeric_limits<double>::min(); 
  double weight_total = 0; 
  vector<double> weights;
  for (size_t i = 0; i < particles.size(); i++)
  {
    weights.push_back(particles[i].weight);
    weight_total += particles[i].weight;
    /* 
    if(weight_max < particles[i].weight)
    {
      weight_max = particles[i].weight;
    }
     */
  }
  
  double weight_max = *max_element(weights.begin(), weights.end());
  
  std::cout << "weight_max " << weight_max << std::endl;
  
  vector<Particle> particles_sampled;
  int N = particles.size();
  
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  
  std::uniform_int_distribution<int>     index_dist(0,N-1);
  std::uniform_real_distribution<double> beta_dist(0.0, weight_max);
  
  int    index = index_dist(gen);
  double beta = 0.0;  
  
  for (int i = 0; i < N; i++)
  {
    beta += beta_dist(gen) * 2.0 ;
    while (beta > weights[index])
    {
      beta -= weights[index];
      index = (index + 1) % N;
    }
    particles_sampled.push_back(particles[index]);
  } 
  particles = particles_sampled;
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