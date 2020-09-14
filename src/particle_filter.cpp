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
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;


using namespace std;



void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  
  
  num_particles = 100;  // TODO: Set the number of particles
      std::default_random_engine generator;
    std::normal_distribution<double> dist1(x, std[0]);
      std::normal_distribution<double> dist2(y, std[1]);
      std::normal_distribution<double> dist3(theta, std[2]);
//   // Resize weights vector based on num_particles
  weights.resize(num_particles);
    
//   // Resize vector of particles
   particles.resize(num_particles);
  
  for(int i=0;i<num_particles;i++){
    Particle tempParticle{0};
  tempParticle.weight=1.0;
    tempParticle.x= dist1(generator);
        tempParticle.y=dist2(generator);
            tempParticle.theta=dist3(generator);
      tempParticle.id=i;
    particles[i]=tempParticle;
    
  }
  is_initialized=true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
      std::default_random_engine generator;
    std::normal_distribution<double> dist1(0, std_pos[0]);
      std::normal_distribution<double> dist2(0, std_pos[1]);
      std::normal_distribution<double> dist3(0, std_pos[2]);
 for(int i=0;i<num_particles;i++){
auto& p=particles[i];
   if (fabs(yaw_rate) < 0.00001) {
 p.x+=velocity*delta_t*cos(p.theta)+dist1(generator);
    p.y+=velocity*delta_t*sin(p.theta)+dist2(generator);
   p.theta+=dist3(generator);
   }
   else{
      p.x+=velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta))+dist1(generator);
    p.y+=velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t))+dist3(generator);
   p.theta+=yaw_rate*delta_t+dist3(generator);
   }
    
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  for(int i=0;i<observations.size();i++){
    int matchedID;
    double minDist{1000000};
    double currDist{-1};
    for(int j=0;j<predicted.size();j++)
    {
      auto pred=predicted[j];
      auto obs=observations[i];
      //currDist=sqrt((obs.x-pred.x)*(obs.x-pred.x)+(obs.y-pred.y)*(obs.y-pred.y));
      currDist=dist(obs.x,obs.y,pred.x,pred.y);

      if(currDist<=minDist){
        minDist=currDist;
        matchedID=pred.id;
      }
    }
    observations[i].id=matchedID;
  }
  

  

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
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
  


       double sum_weights=0;
 for(int i=0;i<num_particles;i++){
  vector<LandmarkObs> predicted;
  vector<LandmarkObs> obsMAP;
   auto& p=particles[i];
   for(int j=0;j<observations.size();j++){
     auto o=observations[j];
        LandmarkObs o_map;
   o_map.x=observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta) + particles[i].x;
   o_map.y= observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta) + particles[i].y;
     //o_map.id=o.id;
     obsMAP.push_back(o_map);
   }
   
	for(int k=0;k<map_landmarks.landmark_list.size();++k){
      LandmarkObs lmrk;
      auto& m=map_landmarks.landmark_list[k];
      double d=dist(m.x_f,m.y_f,p.x,p.y);
      if(d<=sensor_range){
        lmrk.id=m.id_i;
        lmrk.x=m.x_f;
        lmrk.y=m.y_f;
        predicted.push_back(lmrk);
      }
      
    }
   
//  p.x+=velocity*delta_t*cos(p.theta)+dist1(generator);
//     p.y+=velocity*delta_t*sin(p.theta)+dist2(generator);
//    p.theta+=dist3(generator);
  dataAssociation(predicted,  obsMAP);
   p.weight=1.0;
   
      for(int j=0;j<obsMAP.size();j++){
        auto ldmrkID=obsMAP[j].id;
        double pr_x,pr_y;
        for (unsigned int k = 0; k < predicted.size(); k++) {
        	if (predicted[k].id == ldmrkID) {
          		 pr_x = predicted[k].x;
          		 pr_y = predicted[k].y;
        	}
         // std::cout<<pr_x;
        }
        //auto m=map_landmarks.landmark_list[ldmrkID];
         p.weight *= multiv_prob(std_landmark[0], std_landmark[1], obsMAP[j].x, obsMAP[j].y, pr_x,pr_y);
            weights[i] = p.weight;

      }
   //std::cout<<p.weight;
   sum_weights+=p.weight;


  }
     for(int i=0;i<num_particles;i++){
       particles[i].weight=particles[i].weight/sum_weights;
        weights[i] = particles[i].weight;
     }
  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  
  std::default_random_engine generator;
  std::discrete_distribution<int> distribution( weights.begin(), weights.end());
  auto particlesResampled=particles;
  for(int i=0;i<num_particles;i++){
      particlesResampled[i]=particles[distribution(generator)];
  }
particles=particlesResampled;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}