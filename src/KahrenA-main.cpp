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

//********************************************************************
// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s)
{
	auto found_null = s.find("null");
	auto b1 = s.find_first_of("[");
	auto b2 = s.find_first_of("}");
	if (found_null != string::npos) 
	{
		return "";
	} 
	else if (b1 != string::npos && b2 != string::npos) 
	{
	    return s.substr(b1, b2 - b1 + 2);
	}
  	return "";
}

//*********************************************************
double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

//*********************************************************
int ClosestWaypoint(double x, double y, 
					const vector<double> &maps_x, const vector<double> &maps_y)
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

//*********************************************************
int NextWaypoint(double x, double y, double theta, 
					const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

//*******************************************************************
// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
//*******************************************************************
vector<double> getFrenet(double x, double y, double theta, 
						 const vector<double> &maps_x,
						 const vector<double> &maps_y)
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

//*********************************************************
// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, 
					 const vector<double> &maps_x, const vector<double> &maps_y)
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









//*********************************************************
int main() 
{
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
	while (getline(in_map_, line))
	{
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
	
		// Start lane
		int car_lane = 1;

  	// Have a reference velocity to target
		double ref_vel = 4;	// 1.0 mph == 0.623meters/s

	h.onMessage([&ref_vel, &car_lane, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy]		(uWS::WebSocket<uWS::SERVER> ws, 
				char *data, size_t length, uWS::OpCode opCode){

    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') 
		{
      auto s = hasData(data);
      if (s != "") 
	  	{
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") 
				{
          // j[1] is the data JSON object
          
        	// Main car's localization Data
         	double car_x = j[1]["x"];
         	double car_y = j[1]["y"];
					std::cout << "\ncar_x = " << car_x << "\t" << "car_y = " << car_y << "\n";
         	double car_s = j[1]["s"];
         	double car_d = j[1]["d"];
					std::cout << "car_s = " << car_s << "\t" << "car_d = " << car_d << "\n";
         	double car_yaw = j[1]["yaw"];
         	double car_speed = j[1]["speed"];
					std:cout << "car_yaw = " << car_yaw << "\t" << "car_speed = " << car_speed << "\n";

         	// Previous path data given to the Planner
         	auto previous_path_x = j[1]["previous_path_x"];
         	auto previous_path_y = j[1]["previous_path_y"];
			
   				// Previous path's end s and d values 
         	double end_path_s = j[1]["end_path_s"];
         	double end_path_d = j[1]["end_path_d"];
					cout << "end_path_s = " << end_path_s << "\t" << "end_path_d = " << end_path_d << "\n"; 	

         	// Sensor Fusion Data, a list of all other cars on the same side of the road.
					// The data format for each car is: [ id, x, y, vx, vy, s, d]. 
					// The id is a unique identifier for that car. 
					// The x, y values are in global map coordinates, 
					// the vx, vy values are the velocity components, also in reference to the global map. 
					// Finally s and d are the Frenet coordinates for that car
         	vector< vector< double> > sensor_fusion = j[1]["sensor_fusion"];

	//				cout << "sensor_fusion size = " << sensor_fusion.size() << "\n";
	
					double iCar_id, iCar_x, iCar_y, iCar_vx, iCar_vy, iCar_speed, iCar_d, iCar_s;
					double closest_inlane_carId = -1;
					double closest_inlane_car_s;
					double closest_left_carId = -1;
					double closest_left_car_s;
					double closest_right_carId = -1;
					double closest_right_car_s;

					// is this car in my lane? note each lane is 4m wide 
					double lane_min = (double)car_lane * 4;
					double lane_max = (double)(car_lane + 1) * 4;
					cout << "lane_min == "  << lane_min << "\t" << "lane_max == " << lane_max << "\n";

					int slow_down = true;		// assume we can't change lane

					//----------------------------------------
					// Let's find the closest car in my lane 
					//----------------------------------------
					for (int i=0; i<sensor_fusion.size(); i++)
					{
						iCar_id = sensor_fusion[i][0];
						iCar_x = sensor_fusion[i][1];
						iCar_y = sensor_fusion[i][2];
						iCar_vx = sensor_fusion[i][3];
						iCar_vy = sensor_fusion[i][4];
						iCar_speed = sqrt(iCar_vx * iCar_vx + iCar_vy * iCar_vy);
						iCar_s = sensor_fusion[i][5];
						iCar_d = sensor_fusion[i][6]; 

						// If the iCar is not moving ignore or it is behind me
						if( iCar_speed < 1.0 || iCar_s < car_s ) 
						{
//						cout << "ignore carID : " << iCar_id << "\n";
							continue;
						}
//						cout << "\niCar_id = " << iCar_id << "\t" << "ix = " << iCar_x << "\t" << "iy = " <<
//									"\t" << "vx = " << iCar_vx << "\t" << "vy = " << iCar_vy << "speed = " << iCar_speed <<
//									"\t" << "iCar_s = " << iCar_s << "\t" << "iCar_d = " << iCar_d << "\n"; 

						// Is the car within my lane?
						if( (iCar_d > lane_min) && (iCar_d < lane_max)  )
						{
//							cout << "car ID: " << iCar_id << " is in my lane\n";

							// first occurrence of a car in our lane 
							if(closest_inlane_carId == -1)
							{
								closest_inlane_carId = iCar_id;
								closest_inlane_car_s = iCar_s;
							}
							else if	( iCar_s < closest_inlane_car_s) 	// Take the closer one
							{
								closest_inlane_carId = iCar_id;
								closest_inlane_car_s = iCar_s;
							}
						}
						// Is the iCar in my left lane?
						else if ( lane_min >= 4 && (iCar_d > lane_min-4) && (iCar_d < lane_max-4) )
						{
	//						cout << "CarID : " << iCar_id << " is on my left! \n";

							// Determine closest car_to_left 
							if(closest_left_carId == -1) 
							{
								closest_left_carId = iCar_id;
								closest_left_car_s = iCar_s;
							}
							else if (iCar_s < closest_left_car_s)
							{
								closest_left_carId = iCar_id;
								closest_left_car_s = iCar_s;
							}
								
						} // else if
						// Is the iCar in my right lane?
						else if ( lane_min <= 8 && (iCar_d > lane_min+4) && (iCar_d < lane_max+4) )
						{
	//						cout << "CarID : " << iCar_id << " is on my right! \n";

							// Determine closest car to the right
							if(closest_right_carId == -1) 
							{
								closest_right_carId = iCar_id;
								closest_right_car_s = iCar_s;
							}
							else if (iCar_s < closest_right_car_s)
							{
								closest_right_carId = iCar_id;
								closest_right_car_s = iCar_s;
							}
						}
				} // for closest car 

				cout << "closest in lane carID = " << closest_inlane_carId << "\n";
				cout << "closest left lane carID = " << closest_left_carId << "\n";
				cout << "closest right lane carID = " << closest_right_carId << "\n";

				//============================================================		
				//if no cars ahead of us  
				if (closest_inlane_carId == -1)
				{
					cout << "No in-lane cars detected!!\n";
					if (ref_vel < 5)
					{
						cout << "Increasing speed by 1....#1... \n";
						ref_vel += 0.2;			// 0.024m/s == 0.0385 mph :: last value 1.0
					}
					else if( ref_vel < 25)
					{	
						cout << "Increasing speed by 2...#1...\n";
						ref_vel += 2.0;							
					}				
					else if ( ref_vel < 49) 
					{
						cout << "Increasing speed by 1....#2... \n";
						ref_vel += 1.0;			// 0.024m/s == 0.0385 mph 
					}
					else
					{
						cout << "already going @ speed limit \n";
					}
				}			
				else 		// Take the closest car we id'ed
				{
					int i = closest_inlane_carId;	
					iCar_id = sensor_fusion[i][0];
					iCar_x = sensor_fusion[i][1];
					iCar_y = sensor_fusion[i][2];
					iCar_vx = sensor_fusion[i][3];
					iCar_vy = sensor_fusion[i][4];
					iCar_speed = sqrt(iCar_vx * iCar_vx + iCar_vy * iCar_vy);
					iCar_s = sensor_fusion[i][5];
					iCar_d = sensor_fusion[i][6]; 		
	
					// ============================================================
					// given the car is in my lane, is it too close i.e. within 30m?
					if ( iCar_s - car_s > 0 && (iCar_s - car_s < 25) )
					{
						cout << "Car is within 25 meters! and car_lane -1 = " << car_lane-1 << "\n";

						// is there a left lane? 
						if(car_lane -1 >= 0 ) 
						{
							// would I crash into the car if I change to left lane?
							if( closest_left_carId == -1 || abs(sensor_fusion[closest_left_carId][5] - car_s) > 25) 
							{
								cout << "Changing to left lane! \n";
								car_lane -= 1;	// go left
								slow_down = false; 
							}
						}

						// Is there a right lane && we are not changing to left lane?
						if(car_lane +1 < 3 && slow_down == true)
						{
							if( closest_right_carId == -1 || abs(sensor_fusion[closest_right_carId][5] - car_s) > 25) 
							{
								cout << "Changing to right lane! \n";
								car_lane += 1; 	// go right
								slow_down = false; 
							}
						}
						// if we could not change lane for any reason, then slow down
						if(slow_down == true) 
						{
								// can't go left or right... slown down
								cout << "Slowing down ... \n";
								ref_vel -= 1;
						}
					} // there is a car in my lane

					else if(ref_vel < 5) 			// to address the starting jerkiness
					{
						cout << "increasing speed by 1... #3\n";
						ref_vel += 0.2;					// This happens when 
					}				
					else if( ref_vel < 20)   	// no car ahead within range 
					{
						cout << "Increasing speed by 2...#2\n";
						ref_vel += 2;
					}					
					else if ( ref_vel < 49) 
					{
						cout << "Increasing speed by 1...#4 \n";
						ref_vel += 1.0;
					}
					else
					{
						cout << "Already going at speed limit! \n"; 
					}
				}	// else			

				// TODO

//===============================================================================

					int prev_size = previous_path_x.size();
	//				cout << "prev_path_size = " << prev_size << "\n";

					// create a widely spaced (x,y) waypoints, evenly spaced @ 30m
					// later we will interpolate these waypoints with a spline and fill it in with more points 
					// that control speed
					vector <double> ptsx;
					vector <double> ptsy;

					// reference x, y and yaw rate
					// either we will reference the starting point as where the car is or at the previous paths end point
					double ref_x;
					double ref_y;
					double ref_yaw;

					// If previous size is almost empty, use the car as starting reference
					if(prev_size < 2) 
					{
	//					std::cout << "Previous path has: "<< prev_size << " items- simulator has consumed most!!\n";
						ref_x = car_x;
						ref_y = car_y;
						ref_yaw = deg2rad(car_yaw);
	//					cout << "ref_yaw1(deg) = " << ref_yaw << "\n\n";
	
						// Use two points that make the path tangent to the car
						double prev_car_x = car_x - cos(ref_yaw); 
						double prev_car_y = car_y - sin(ref_yaw);
	//					std::cout << "car_x = " << car_x << "\t" << "prev_car_x = " << prev_car_x << "\n";
	//					std::cout << "car_y = " << car_y << "\t" << "prev_car_y = " << prev_car_y << "\n";
	
						ptsx.push_back(prev_car_x);
						ptsx.push_back(car_x);

						ptsy.push_back(prev_car_y);
						ptsy.push_back(car_y);
			
					}
					else   // use the previous path's end point as starting reference
					{
	//					std::cout << "Remaining items in previous path = " << prev_size << "\n";

						// Redefine reference state as previous path end point
						ref_x = previous_path_x[prev_size-1];
						ref_y = previous_path_y[prev_size-1];
				
						double ref_x_prev = previous_path_x[prev_size-2];
						double ref_y_prev = previous_path_y[prev_size-2];

						ref_yaw = atan2( ref_y-ref_y_prev, ref_x-ref_x_prev);		

						// Use two points that make the path tangent to the previous path's end point
						ptsx.push_back(ref_x_prev);
						ptsx.push_back(ref_x);
						ptsy.push_back(ref_y_prev);
						ptsy.push_back(ref_y);
						
					}

					//----------------------------------------------------------------------------
					// Create /add evenly 30m spaced points ahead of the starting reference
					vector<double> next_wp0 = getXY(car_s + 60, (2+4*car_lane), map_waypoints_s, 
																						map_waypoints_x, map_waypoints_y);
					vector<double> next_wp1 = getXY(car_s + 90, (2+4*car_lane), map_waypoints_s, 
																						map_waypoints_x, map_waypoints_y);
					vector<double> next_wp2 = getXY(car_s + 120, (2+4*car_lane), map_waypoints_s, 
																						map_waypoints_x, map_waypoints_y);

					ptsx.push_back(next_wp0[0]);
					ptsx.push_back(next_wp1[0]);
					ptsx.push_back(next_wp2[0]);
	//				std::cout << "ptsx.pushback = " << next_wp0[0] << " " << next_wp1[0] << " " << next_wp2[0] << "\n";

					ptsy.push_back(next_wp0[1]);
					ptsy.push_back(next_wp1[1]);
					ptsy.push_back(next_wp2[1]);
	//				std::cout << "ptsy.pushback = " << next_wp0[1] << " " << next_wp1[1] << " " << next_wp2[1] << "\n";

					//----------------------------------------------------------------------------
					// shift car reference angle to 0 degrees. basically switch to car-coordinates
					//----------------------------------------------------------------------------
					for(int i=0; i< ptsx.size(); i++)
					{
						double delta_x = ptsx[i] - ref_x;
						double delta_y = ptsy[i] - ref_y;

						ptsx[i] = (delta_x * cos(0-ref_yaw)) - (delta_y * sin(0-ref_yaw));
						ptsy[i] = (delta_x * sin(0-ref_yaw)) + (delta_y * cos(0-ref_yaw));
	//					cout << " Car coord ptsx[i] = " << ptsx[i] << "\t" << "ptsy[i] = " << ptsy[i] << "\n";
					}

					// create spline
					tk::spline s;
	
					//set (x,y) points to spline 
					s.set_points(ptsx, ptsy);
	
					//-------------------------------------------------------------------
         	vector<double> next_x_vals;
         	vector<double> next_y_vals;

					int start_of_prev_path = 0;

					// If we are slowing down, we need to take fewer points instead of all
					if(slow_down == true && prev_size > 2) 
					{
						start_of_prev_path = prev_size - 2; 
					}

					// Start with ALL of the previous path points from last time
					for (int i=0; i < prev_size; i++)
					{
						next_x_vals.push_back(previous_path_x[i]);
						next_y_vals.push_back(previous_path_y[i]);
					}
			
					// Calculate how to break up spline points so that we travel at our desired reference velocity
					double target_x = 30.0;
					double target_y = s(target_x);
					double target_dist = sqrt ((target_x)*(target_x) + (target_y)*(target_y) );
	//				std::cout << "target_x = " << target_x << "\t" << "target_y = " << target_y << "\t" 
	//														  << "target_dist = " << target_dist << "\n";

					double x_add_on = 0;

					// Fill up the rest of the path planner after filling it with previous points. 
					// Here we will always put out 50 points
					for(int i = 0; i < 50-prev_size; i++)
					{
						double N = (target_dist / (0.02 * ref_vel / 2.24) );   // check out the hand-drawn diagram
//					std::cout << "N = " << N << "\n";

						double x_point = x_add_on + (target_x)/N;				// 2.24 for mph to meterps
						double y_point = s(x_point);
//					std::cout << "x_point = " << x_point << "\t" << "y_point = " << y_point << "\n";			

						x_add_on = x_point;
//						cout << "x_add_on = " << x_add_on << "\n";

						double x_ref = x_point;
						double y_ref = y_point;

						//---------------------------------------------------------
						// rotate back to World Map coordinates 
						//---------------------------------------------------------
						x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw) );
						y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw) );
//						std::cout << "x_point = " << x_point << "\t" << "y_point = " << y_point << "\n";			

						x_point += ref_x;
						y_point += ref_y;
	//					std::cout << "next_xy_vals: x_point = " << x_point << "\t" << "y_point = " << y_point << "\n";			

						next_x_vals.push_back(x_point);
						next_y_vals.push_back(y_point);

					}	// for

//					cout << "\n\n"; 		

					// end of TODO 

    		  json msgJson;
					msgJson["next_x"] = next_x_vals;
    		  msgJson["next_y"] = next_y_vals;

    		  auto msg = "42[\"control\","+ msgJson.dump()+"]";

    		  //this_thread::sleep_for(chrono::milliseconds(1000));
    		  ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

    		} // if event is telemetry
    	} // s!=0
			else
			{
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
    	}
  	} // msg==42
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
