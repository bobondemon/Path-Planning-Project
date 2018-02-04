#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

/*
*****
***** MY Helper function
*****
*/
const double LARGEDOUBLE = 99999.0;
const double SPEEDLIMIT = 47.5;
// Input d, return which lane index, 0, 1, 2
int getLaneIdx(const double d)
{
  int laneIdx = int(d/4.0);
  if (d<0 || laneIdx>2)
    laneIdx = -1;
  return laneIdx;
}
// As name refers
double EucledianDist(const double x1, const double y1, const double x2, const double y2)
{
  return sqrt(pow(fabs(x1-x2),2) + pow(fabs(y1-y2),2));
}
// Calculate the smallest distance between two trajectories
double EucledianDistBetweenTrajectories(const vector<double> &x1, const vector<double> &y1, const vector<double> &x2, const vector<double> &y2)
{
  double min_dist = LARGEDOUBLE;
  for (int i=0; i<x1.size(); i++)
  {
    double dist = EucledianDist(x1[i],y1[i],x2[i],y2[i]);
    if (dist<min_dist)
      min_dist = dist;
  }
  return min_dist;
}
// Check collision for tow trajectories
bool isCollision(const vector<double> &x1, const vector<double> &y1, const vector<double> &x2, const vector<double> &y2)
{
  const double car_radius = 2.0;
  const double additional_buffer = 3.0;  // Make it safer to do lane changing
  for (int i=0; i<x1.size(); i++)
  {
    double dist = EucledianDist(x1[i],y1[i],x2[i],y2[i]);
    if (dist < car_radius*2+additional_buffer)
      return true;
  }
  return false;
}
// Return true if two trajectories has intersections within radius buffer_radius
bool isNear(const vector<double> &x1, const vector<double> &y1, const vector<double> &x2, const vector<double> &y2)
{
  const double buffer_radius = 30.0;
  for (int i=0; i<x1.size(); i++)
  {
    double dist = EucledianDist(x1[i],y1[i],x2[i],y2[i]);
    if (dist < buffer_radius*2)
      return true;
  }
  return false;
}
// Return the car index that is closest (ahead/behind) of ego car at check_lane
// E.g.: Assume ego car is at lane 1, and we are going to LCR, then check_lane=2.
//       position='ahead' will give us the closest car index that is ahead of ego car at lane 2.
//       If no ahead car at lane 2, then it'll return -1
int getClosestCarIdx(const std::vector<double> &car_data, const std::vector<vector<double>> &sensor_fusion, const int check_lane, const string position)
{
  int rtnIdx = -1;
  double minDist = LARGEDOUBLE;
  double car_x = car_data[0];
  double car_y = car_data[1];
  double car_s = car_data[2];
  for (int i=0; i<sensor_fusion.size(); i++)
  {
    float check_car_d = sensor_fusion[i][6];
    int check_car_lane = getLaneIdx(check_car_d);
    if (check_car_lane != check_lane)
      continue;

    double check_car_s = sensor_fusion[i][5];
    double dist = fabs(check_car_s - car_s);
    if (dist < minDist)
    {
      if (position.compare("behind") == 0 && check_car_s < car_s)
      {
        rtnIdx = i;
        minDist = dist;
      }
      if (position.compare("ahead") == 0 && check_car_s > car_s)
      {
        rtnIdx = i;
        minDist = dist;
      }
    }
  }
  return rtnIdx;
}
// Input the lane index of ego car, then return the possible actions
vector<string> successor_states_IMPL(const int ego_lane)
{
  // Provides the possible next states given the current state for the FSM
  int lanes_available = 3;  // in this simulation we only has three lanes
  vector<string> states;
  states.push_back("KL");
  if (ego_lane != 0)
    states.push_back("LCL");
  if (ego_lane != lanes_available - 1)
    states.push_back("LCR");
  return states;
}
// Cost functions of speed efficiency. A larger speed produces a smaller cost.
double speed_efficiency_cost(const double v)
{
  // Simply a linear cost from 1 to 0, where speed is from 0 to SPEEDLIMIT
  return 1.0-v/SPEEDLIMIT;
}
// Penalize close distance of two trajectories (for safety issue)
double how_close_cost(const vector<double> &x1, const vector<double> &y1, const vector<double> &x2, const vector<double> &y2)
{
  double rtnCost;
  double min_dist = LARGEDOUBLE;
  for (int i=0;i<x1.size();i++)
  {
    double dist = EucledianDist(x1[i],y1[i],x2[i],y2[i]);
    if (dist < min_dist)
      min_dist = dist;
  }
  // distance less than 4.0 will has cost 1
  // otherwirse, cost = (4.0/dist)^2. dist=12 will has cost 0.11111
  // otherwirse, cost = (4.0/dist)^1.1. dist=12 will has cost 0.298
  if (min_dist<=4.0)
    rtnCost = 1.0;
  else
    rtnCost = pow((4.0/min_dist),2);
  return rtnCost;
}
// Predict a car's trajectory with assuming constant velocity
void simplePrediction(const std::vector<vector<double>> &sensor_fusion, const int idx, const int num_of_pred,
  vector<double> &rtn_next_x_vals, vector<double> &rtn_next_y_vals)
{
  double curx = sensor_fusion[idx][1];
  double cury = sensor_fusion[idx][2];
  rtn_next_x_vals.push_back(curx);
  rtn_next_y_vals.push_back(cury);
  double car_vx = sensor_fusion[idx][3];
  double car_vy = sensor_fusion[idx][4];
  for (int i=1;i<num_of_pred;i++)
  {
    curx += .02*car_vx;
    rtn_next_x_vals.push_back(curx);

    cury += .02*car_vy;
    rtn_next_y_vals.push_back(cury);
  }
}

/*
*****
***** Trajectory Generation Functions
*****
*/
// Important arguments are: "sd_points" and "ref_vel", then this function will generate the tragjectory based on target s, d and velocity
void genTrajectory(const vector<double> &map_waypoints_x, const vector<double> &map_waypoints_y,const vector<double> &map_waypoints_s,
  const std::vector<double> &in_previous_path_x, const std::vector<double> &in_previous_path_y, const std::vector<double> &car_data,
  const std::vector<vector<double>> &sd_points, const double ref_vel, vector<double> &next_x_vals, vector<double> &next_y_vals)
{
  const int trajectory_len = 75;  // the trajectory point number that we are going to generate
  const int num_reuse_prev_traj_pts = 30; // number of points that we are going to use in previous trajectory
  vector<double> previous_path_x;
  vector<double> previous_path_y;
  // Preserve num_reuse_prev_traj_pts points of previous trajectories
  int original_prev_size = in_previous_path_x.size();
  int prev_size = original_prev_size;
  if (original_prev_size > num_reuse_prev_traj_pts)
    prev_size = num_reuse_prev_traj_pts;
  for (int i=0; i<prev_size; i++)
  {
    previous_path_x.push_back(in_previous_path_x[i]);
    previous_path_y.push_back(in_previous_path_y[i]);
  }
  // cout << original_prev_size << ",  " << previous_path_y.size() << endl;
  double car_x = car_data[0];
  double car_y = car_data[1];
  double car_s = car_data[2];
  double car_d = car_data[3];
  double car_yaw = car_data[4];
  double car_speed = car_data[5];
  // ===== Using spline to find out the smoothed path
  vector<double> ptsx;
  vector<double> ptsy;
  // reference x, y, yaw states
  // either we'll reference the starting point as where the car is or at the previous paths end point
  double ref_x = car_x;
  double ref_y = car_y;
  double ref_yaw = deg2rad(car_yaw);
  // if previous size is almost empty, use the car as starting reference
  if (prev_size < 2)
  {
    // Use two points that make the path tangent to the car
    double prev_car_x = car_x - cos(ref_yaw);
    double prev_car_y = car_y - sin(ref_yaw);

    ptsx.push_back(prev_car_x);
    ptsx.push_back(car_x);

    ptsy.push_back(prev_car_y);
    ptsy.push_back(car_y);
  }
  // use the previous path's end point as starting reference
  else
  {
    // redefine reference state as previous path end point
    ref_x = previous_path_x[prev_size-1];
    ref_y = previous_path_y[prev_size-1];

    double prev_ref_x = previous_path_x[prev_size-2];
    double prev_ref_y = previous_path_y[prev_size-2];
    ref_yaw = atan2(ref_y - prev_ref_y, ref_x - prev_ref_x);

    // Use two points that make the path tangent to the previous path's end point
    ptsx.push_back(prev_ref_x);
    ptsx.push_back(ref_x);

    ptsy.push_back(prev_ref_y);
    ptsy.push_back(ref_y);
  }
  // So far we push two x and two y in pts
  // In Frenet add evenly 30m second points ahead of the starting reference
  for (int i=0;i<sd_points.size();i++)
  {
    vector<double> next_wp = getXY(car_s+sd_points[i][0],sd_points[i][1],map_waypoints_s,map_waypoints_x,map_waypoints_y);
    ptsx.push_back(next_wp[0]);
    ptsy.push_back(next_wp[1]);
  }

  // Do global to local coordinates transformation
  for (int i=0; i<ptsx.size(); i++)
  {
    // shift car reference coordinates to (0, 0)
    double shift_x = ptsx[i]-ref_x;
    double shift_y = ptsy[i]-ref_y;
    ptsx[i] = (shift_x*cos(0-ref_yaw) - shift_y*sin(0-ref_yaw));
    ptsy[i] = (shift_x*sin(0-ref_yaw) + shift_y*cos(0-ref_yaw));
  }

  // create a spline
  tk::spline s;
  // set (x,y) points to the spline
  s.set_points(ptsx,ptsy);

  // Start with all of the previous path points from last time
  for (int i=0; i<prev_size; i++)
  {
    next_x_vals.push_back(previous_path_x[i]);
    next_y_vals.push_back(previous_path_y[i]);
  }

  // Calculate how to break up spline points so that we travel at our desired referene velocity
  double target_x = 30.0;
  double target_y = s(target_x);
  double target_dist = sqrt(target_x*target_x + target_y*target_y);
  double x_add_on = 0;

  // Fill up the rest of our path planner after filling it with previous points, here we'll always output 50 points
  for (int i=1; i<=trajectory_len-prev_size; i++)
  {
    double N = target_dist/(.02*ref_vel/2.24); // 2.24 is because ref_vel is in mph, and map coordinates are in meter
    double x_point = x_add_on+(target_x)/N;
    double y_point = s(x_point);

    x_add_on = x_point;

    double x_ref = x_point;
    double y_ref = y_point;
    // rotate back to map coordinate
    x_point = x_ref*cos(ref_yaw) - y_ref*sin(ref_yaw);
    y_point = x_ref*sin(ref_yaw) + y_ref*cos(ref_yaw);

    x_point += ref_x;
    y_point += ref_y;

    next_x_vals.push_back(x_point);
    next_y_vals.push_back(y_point);
  }
}

// Return trajectory that keep in current lane. Slow down if too close to the car infront of us, otherwise, speed up to target speed
double keep_lane_trajectory(const vector<double> &map_waypoints_x, const vector<double> &map_waypoints_y,const vector<double> &map_waypoints_s,
  const std::vector<double> &previous_path_x, const std::vector<double> &previous_path_y,
  const string state, const std::vector<double> &car_data, const std::vector<vector<double>> &sensor_fusion,
  double *ref_vel_ptr, vector<double> &rtn_next_x_vals, vector<double> &rtn_next_y_vals)
{
  const int check_car_trajectory_len = 75;  // the trajectory point number that we are going to generate for OTHER CARS
  double rtnCost = LARGEDOUBLE;
  double car_x = car_data[0];
  double car_y = car_data[1];
  double car_s = car_data[2];
  double car_d = car_data[3];
  double car_yaw = car_data[4];
  double car_speed = car_data[5];
  double end_path_s = car_data[6];
  double end_path_d = car_data[7];
  
  int prev_size = previous_path_x.size();
  bool too_close = false;
  int ego_lane = getLaneIdx(car_d);

  double check_car_speed = SPEEDLIMIT;
  double nearest_s_ahead = LARGEDOUBLE;
  // ===== stay away from the car in front of us
  for (int i=0; i<sensor_fusion.size(); i++)
  {
    // car is in my lane
    float d = sensor_fusion[i][6];
    int check_car_lane = getLaneIdx(d);
    if (check_car_lane == ego_lane)
    {
      double vx = sensor_fusion[i][3];
      double vy = sensor_fusion[i][4];
      double check_car_s = sensor_fusion[i][5];
      double this_car_speed  = sqrt(vx*vx+vy*vy);
      if (check_car_s<nearest_s_ahead)
      {
        nearest_s_ahead = check_car_s;
        check_car_speed  = this_car_speed;
      }
      // Check if we have a collision in the future.
      check_car_s += (prev_size*.02*this_car_speed);
      // check s values greater than mine and s gap
      if ( (check_car_s>end_path_s) && (check_car_s-end_path_s<30) )
      {
        too_close = true;
        // must break, otherwise the check_car_speed may get wrong value
        // we need check_car_speed be the speed of the "too_close" car infront of us
        // and use check_car_speed to calculate the speed_efficiency_cost
        break;
      }
    }
  }

  if (too_close)
  {
    *ref_vel_ptr -= .2;  // .224
  }
  else if (*ref_vel_ptr < SPEEDLIMIT)
  {
    *ref_vel_ptr += .2;  // .224
  }
  // Generate sd_points for KL trajectory
  std::vector<vector<double>> sd_points;
  for (int i=1;i<=3;i++)
  {
    std::vector<double> sd = {i*30.0, (2.0+4.0*ego_lane)};
    sd_points.push_back(sd);
  }

  // Call tragjectory generation
  genTrajectory(map_waypoints_x, map_waypoints_y, map_waypoints_s, previous_path_x, previous_path_y,
              car_data, sd_points, *ref_vel_ptr, rtn_next_x_vals, rtn_next_y_vals);

  // Predict the trajectories of the car in front of the ego lane
  int carIdxAhead = getClosestCarIdx(car_data, sensor_fusion, ego_lane, "ahead");
  std::vector<double> aheadCarTrajectory_x;
  std::vector<double> aheadCarTrajectory_y;
  if (carIdxAhead!=-1) // Have Car Ahead
  {
    simplePrediction(sensor_fusion, carIdxAhead, check_car_trajectory_len, aheadCarTrajectory_x, aheadCarTrajectory_y);
  }
  // Calculate the cost
  double speed_cost = speed_efficiency_cost(check_car_speed);
  if (!isNear(aheadCarTrajectory_x, aheadCarTrajectory_y, rtn_next_x_vals, rtn_next_y_vals))
    speed_cost = 0;
  double change_lane_cost = 0;
  rtnCost = change_lane_cost + speed_cost;
  return rtnCost;
}

// Return trajectory that can safely change lane, otherwise, return empty trajectory
double lane_change_trajectory(const vector<double> &map_waypoints_x, const vector<double> &map_waypoints_y,const vector<double> &map_waypoints_s,
  const std::vector<double> &previous_path_x, const std::vector<double> &previous_path_y,
  const string state, const std::vector<double> &car_data, const std::vector<vector<double>> &sensor_fusion,
  double *ref_vel_ptr, vector<double> &rtn_next_x_vals, vector<double> &rtn_next_y_vals)
{
  double rtnCost = LARGEDOUBLE;
  const int check_car_trajectory_len = 75;  // the trajectory point number that we are going to generate for OTHER CARS

  double car_x = car_data[0];
  double car_y = car_data[1];
  double car_s = car_data[2];
  double car_d = car_data[3];
  double car_yaw = car_data[4];
  double car_speed = car_data[5];
  double end_path_s = car_data[6];
  double end_path_d = car_data[7];

  int ego_lane = getLaneIdx(car_d);
  int intend_lane = ego_lane;
  if (state.compare("LCL")==0)
    intend_lane --;
  else if (state.compare("LCR")==0)
    intend_lane++;
  
  // ===== stay away from the car in front of us
  int prev_size = previous_path_x.size();
  bool too_close = false;

  double ego_lane_ahead_car_speed = SPEEDLIMIT;  // this is the speed of car which is in front of ego car at ego lane
  double nearest_s_ahead = LARGEDOUBLE;
  for (int i=0; i<sensor_fusion.size(); i++)
  {
    // car is in my lane
    float d = sensor_fusion[i][6];
    int check_car_lane = getLaneIdx(d);
    if (check_car_lane == ego_lane)
    {
      double vx = sensor_fusion[i][3];
      double vy = sensor_fusion[i][4];
      double check_car_s = sensor_fusion[i][5];
      double this_car_speed  = sqrt(vx*vx+vy*vy);
      if (check_car_s<nearest_s_ahead)
      {
        nearest_s_ahead = check_car_s;
        ego_lane_ahead_car_speed  = this_car_speed;
      }
      // Check if we have a collision in the future.
      check_car_s += (prev_size*.02*this_car_speed);
      // check s values greater than mine and s gap
      if ( (check_car_s>end_path_s) && (check_car_s-end_path_s<30) )
      {
        too_close = true;
        break;
      }
    }
  }

  // Predict the trajectories of the cars ahead/behind at the intended lane
  int carIdxBehind = getClosestCarIdx(car_data, sensor_fusion, intend_lane, "behind");
  std::vector<double> behindCarTrajectory_x;
  std::vector<double> behindCarTrajectory_y;
  double behind_car_speed  = 0;
  if (carIdxBehind!=-1) // Have Car Behind
  {
    simplePrediction(sensor_fusion, carIdxBehind, check_car_trajectory_len, behindCarTrajectory_x, behindCarTrajectory_y);
    double behind_car_vx = sensor_fusion[carIdxBehind][3];
    double behind_car_vy = sensor_fusion[carIdxBehind][4];
    behind_car_speed  = sqrt(behind_car_vx*behind_car_vx+behind_car_vy*behind_car_vy);
  }
  int carIdxAhead = getClosestCarIdx(car_data, sensor_fusion, intend_lane, "ahead");
  std::vector<double> aheadCarTrajectory_x;
  std::vector<double> aheadCarTrajectory_y;
  double ahead_car_speed  = SPEEDLIMIT;
  bool ahead_car_close = false;
  if (carIdxAhead!=-1) // Have Car Ahead
  {
    simplePrediction(sensor_fusion, carIdxAhead, check_car_trajectory_len, aheadCarTrajectory_x, aheadCarTrajectory_y);
    double ahead_car_vx = sensor_fusion[carIdxAhead][3];
    double ahead_car_vy = sensor_fusion[carIdxAhead][4];
    ahead_car_speed  = sqrt(ahead_car_vx*ahead_car_vx+ahead_car_vy*ahead_car_vy);
  }
  // Change velocity to fit the intended lane, using ego_lane_ahead_car_speed, ahead_car_speed
  // If ahead_car at intended lane is far away, then we can speed up to change lane
  if (carIdxAhead!=-1 || !isNear(aheadCarTrajectory_x, aheadCarTrajectory_y, rtn_next_x_vals, rtn_next_y_vals))
  {
    if (*ref_vel_ptr < SPEEDLIMIT)
    {
      if (too_close)  // we have a close car in font of us at our lane, so not speed up too much
        *ref_vel_ptr += .1;
      else
        *ref_vel_ptr += .2;  // .224
    }
  }

  // Generate sd_points for intended lane
  std::vector<vector<double>> sd_points;
  double intend_d = (2.0+4.0*intend_lane);
  std::vector<double> dvec = {intend_d,intend_d,intend_d};
  // std::vector<double> dvec = {car_d+3.0*(intend_d-car_d)/4.0,intend_d,intend_d};
  for (int i=1;i<=3;i++)
  {
    std::vector<double> sd = {i*50.0, dvec[i-1]};
    sd_points.push_back(sd);
  }
  // Call tragjectory generation
  genTrajectory(map_waypoints_x, map_waypoints_y, map_waypoints_s, previous_path_x, previous_path_y,
              car_data, sd_points, *ref_vel_ptr, rtn_next_x_vals, rtn_next_y_vals);
  
  // Check Collisions
  if (carIdxBehind!=-1) // Have Car Behind
  {
    if (isCollision(behindCarTrajectory_x, behindCarTrajectory_y, rtn_next_x_vals, rtn_next_y_vals))
    {
      cout << "***** Predicting Collision with BEHIND car at " << (((ego_lane-intend_lane)<0)?"Right":"Left") << " lane\n";
      rtn_next_x_vals.clear();  // having collision, then return empty trajectory
      rtn_next_y_vals.clear();  // having collision, then return empty trajectory
    }
  }

  if (carIdxAhead!=-1) // Have Car Ahead
  {
    // Check if we have a ahead car in horizon at intended lane
    if (isNear(aheadCarTrajectory_x, aheadCarTrajectory_y, rtn_next_x_vals, rtn_next_y_vals))
      ahead_car_close = true;
    // Check collision of this trajectory
    if (isCollision(aheadCarTrajectory_x, aheadCarTrajectory_y, rtn_next_x_vals, rtn_next_y_vals))
    {
      cout << "***** Predicting Collision with AHEAD car at " << (((ego_lane-intend_lane)<0)?"Right":"Left") << " lane\n";
      rtn_next_x_vals.clear();  // having collision, then return empty trajectory
      rtn_next_y_vals.clear();  // having collision, then return empty trajectory
    }
  }

  // Calculate the cost
  double speed_cost = speed_efficiency_cost(ahead_car_speed);
  double dist_factor = 1.0;
  if (carIdxAhead!=-1)  // Descrease the speed_cost by a factor if having ahead car at intended lane
  {
    double dist = EucledianDistBetweenTrajectories(aheadCarTrajectory_x, aheadCarTrajectory_y, rtn_next_x_vals, rtn_next_y_vals);
    dist_factor = 1.0 - 1.0/(1.0+exp((60.0-dist)*1.0));
  }
  // cout << "speed_cost=" << speed_cost;
  // speed_cost *= dist_factor;
  // cout << "; dist_factor=" << dist_factor << "; speed_cost*facotr=" << speed_cost << endl;
  
  // if (!ahead_car_close)
  //   speed_cost = 0;

  double change_lane_cost = 0.10;
  double how_close_behind_car_cost = 0.0;
  if (behindCarTrajectory_x.size()!=0)
    how_close_behind_car_cost = how_close_cost(behindCarTrajectory_x,behindCarTrajectory_y,rtn_next_x_vals,rtn_next_y_vals);
  double how_close_ahead_car_cost = 0.0;
  if (aheadCarTrajectory_x.size()!=0)
    how_close_ahead_car_cost = how_close_cost(aheadCarTrajectory_x,aheadCarTrajectory_y,rtn_next_x_vals,rtn_next_y_vals);
  rtnCost = change_lane_cost + speed_cost;
  rtnCost += (how_close_behind_car_cost<how_close_ahead_car_cost)?how_close_ahead_car_cost:how_close_behind_car_cost;  // chose max to add up
  return rtnCost;
}

// Generate State Trajectory
// Given a possible state, GST will generate the appropriate trajectory to realize that state action.
double GST_IMPL(const vector<double> &map_waypoints_x, const vector<double> &map_waypoints_y,const vector<double> &map_waypoints_s,
  const std::vector<double> &in_previous_path_x, const std::vector<double> &in_previous_path_y,
  const string state, const std::vector<double> &car_data, const std::vector<vector<double>> &sensor_fusion,
  double *ref_vel_ptr, vector<double> &rtn_next_x_vals, vector<double> &rtn_next_y_vals)
{
  double rtnCost = LARGEDOUBLE;

  if (state.compare("KL") == 0) {
    rtnCost = keep_lane_trajectory(map_waypoints_x, map_waypoints_y, map_waypoints_s,
      in_previous_path_x, in_previous_path_y,
      state, car_data, sensor_fusion,
      ref_vel_ptr, rtn_next_x_vals, rtn_next_y_vals);
  } else if (state.compare("LCL") == 0 || state.compare("LCR") == 0) {
    rtnCost = lane_change_trajectory(map_waypoints_x, map_waypoints_y, map_waypoints_s,
      in_previous_path_x, in_previous_path_y,
      state, car_data, sensor_fusion,
      ref_vel_ptr, rtn_next_x_vals, rtn_next_y_vals);
  }
  return rtnCost;
}

/*
*****
***** Main function
*****
*/
int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  // Have a reference velocity to target
  double ref_vel = 0;  // mph
  // double ref_vel = 49.5;  // mph
  string prev_state = "KL";

  // ===== Define the function pointers, (Generate State Trajectory)
  double (*GST)(const vector<double> &map_waypoints_x, const vector<double> &map_waypoints_y,const vector<double> &map_waypoints_s,
    const std::vector<double> &in_previous_path_x, const std::vector<double> &in_previous_path_y,
    const string state, const std::vector<double> &car_data, const std::vector<vector<double>> &sensor_fusion,
    double *ref_vel_ptr, vector<double> &next_x_vals, vector<double> &next_y_vals);
  GST = GST_IMPL;

  vector<string> (*successor_states)(const int lane_idx);
  successor_states = successor_states_IMPL;

  h.onMessage([&prev_state, &successor_states, &GST, &ref_vel, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // cout << pfun(1.5,2.5) << endl;
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

            int ego_lane = getLaneIdx(car_d);

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
            std::vector<double> car_data={car_x,car_y,car_s,car_d,car_yaw,car_speed,end_path_s,end_path_d};

            // ===== Behavior Planning with Finite State Machine
            // Choose the states with the lowest cost based on the predictions of other car's info in sensor_fusion
            vector<string> possible_states = successor_states(ego_lane);
            std::vector<vector<double>> trajectories_x;
            std::vector<vector<double>> trajectories_y;
            std::vector<string> trajectories_states;
            std::vector<double> costs;
            std::vector<double> ref_vels;

            for (int i=0; i<possible_states.size(); i++)
            {
              vector<double> trajectory_x;
              vector<double> trajectory_y;
              double state_ref_vel = ref_vel;
              // Generate State Trajectory (GST)
              double cost = GST(map_waypoints_x, map_waypoints_y, map_waypoints_s, previous_path_x, previous_path_y,
                possible_states[i], car_data, sensor_fusion,
                &state_ref_vel, trajectory_x, trajectory_y);
              if (trajectory_x.size()!=0)  // In GST, if collisions predicted, then trajectory_x/y will be empty
              {
                trajectories_x.push_back(trajectory_x);
                trajectories_y.push_back(trajectory_y);
                trajectories_states.push_back(possible_states[i]);
                ref_vels.push_back(state_ref_vel);
                // Add state transition cost, prefer keeping at same state
                double state_transition_cost = 0.0;
                if (prev_state.compare(possible_states[i])!=0)
                {
                  state_transition_cost = 0.05;
                  if ( (prev_state.compare("LCL")==0) && (possible_states[i].compare("LCR")==0) )
                    state_transition_cost = 0.2;
                  if ( (prev_state.compare("LCR")==0) && (possible_states[i].compare("LCL")==0) )
                    state_transition_cost = 0.2;
                }
                cost += state_transition_cost;
                costs.push_back(cost);
              }
            }

            // Debug messages
            for (int i=0;i<costs.size();i++)
            {
              cout << possible_states[i] << "=" << costs[i] << "; ";
            }
            cout << endl;

            // Choose the lowest costs
            int min_idx = 0;
            double min_cost = costs[min_idx];
            for (int i=1;i<costs.size();i++)
            {
              if (costs[i]<min_cost)
              {
                min_idx = i;
                min_cost = costs[i];
              }
            }

            cout << "Previous state " << prev_state << "; Take action " << trajectories_states[min_idx] << endl;
            prev_state = trajectories_states[min_idx];

            ref_vel = ref_vels[min_idx];
            vector<double> next_x_vals = trajectories_x[min_idx];
            vector<double> next_y_vals = trajectories_y[min_idx];

            json msgJson;

          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
