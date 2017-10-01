#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

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
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

Eigen::VectorXd globalKinematic(Eigen::VectorXd state,
                                Eigen::VectorXd actuators, double dt) {
  Eigen::VectorXd next_state(state.size());

  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double delta = actuators[0];
  double a = actuators[1];
  double Lf = 2.67;

  next_state[0] = x + v*cos(psi)*dt;
  next_state[1] = y + v*sin(psi)*dt;
  next_state[2] = psi + v/Lf * delta * dt;
  next_state[3] = v + a*dt;

  return next_state;
}



int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];

          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */
          // reference to the quiz in the class

//          //---- consider the latency -----
          Eigen::VectorXd state_temp(4);
          Eigen::VectorXd actuators_temp(2);
          state_temp << px, py, psi, v;
          actuators_temp << 0, 0;
          auto next_state_latency = globalKinematic(state_temp, actuators_temp, 0.1);

          Eigen::VectorXd ptsx_eigen(ptsx.size());
          Eigen::VectorXd ptsy_eigen(ptsy.size());
          px = next_state_latency[0];
          py = next_state_latency[1];
          //---- Transformed waypoints into the vehicle frame -----
          for(unsigned int k=0; k< ptsx.size(); k++){
            double x = ptsx[k] - px;
            double y = ptsy[k] - py;

            ptsx_eigen[k] = cos(psi)*x + sin(psi)*y;
            ptsy_eigen[k] = -sin(psi)*x + cos(psi)*y;
//              ptsx_eigen[k] = ptsx[k];
//              ptsy_eigen[k] = ptsy[k];

          }

//          for(unsigned int xi=0; xi<ptsx.size(); xi++){
//            ptsx_eigen[xi] = ptsx[xi];
//          }
//          for(unsigned int yi=0; yi<ptsx.size(); yi++){
//            ptsy_eigen[yi] = ptsy[yi];
//          }
          auto coeffs = polyfit(ptsx_eigen, ptsy_eigen, 3);

          // The cross track error is calculated by evaluating at polynomial at x, f(x)
          // and subtracting y.
          double predicted_y = polyeval(coeffs, ptsx_eigen[0]); //px?
          double cte = predicted_y - ptsy_eigen[0]; //py?
          // Due to the sign starting at 0, the orientation error is -f'(x).
          // derivative of coeffs[0] + coeffs[1] * x -> coeffs[1]
          double epsi = psi - atan(coeffs[1]);

          Eigen::VectorXd state(6);
//          state << px, py, psi, v, cte, epsi;
          // [x, y, psi, v]

          state << 0,0, 0, v, cte, epsi;
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          double steer_value=0;
          double throttle_value=0;

          auto vars = mpc.Solve(state, coeffs);

          for (unsigned int n=0; n<(unsigned int)(vars.size()/8); n++){
            mpc_x_vals.push_back(vars[0+8*n]);
            mpc_y_vals.push_back(vars[1+8*n]);
          }
          steer_value = -vars[6];
          throttle_value = vars[7];

//          for (unsigned int k=0; k<1;k++){
//            auto vars = mpc.Solve(state, coeffs);
//            mpc_x_vals.push_back(vars[0]);
//            mpc_y_vals.push_back(vars[1]);
//            state << ptsx_eigen[k],ptsy_eigen[k], vars[2], vars[3], vars[4], vars[5];
//            if(k==0){
//              steer_value = -vars[6];
//              throttle_value = vars[7];
//            }
//
//          }



          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].

          steer_value = steer_value/deg2rad(25);

          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;


          //Display the MPC predicted trajectory 


          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          for(unsigned int k=0; k< 20; k++){
            predicted_y = polyeval(coeffs, k*2.0);
            next_x_vals.push_back(k*2.0);
            next_y_vals.push_back(predicted_y);
          }



          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
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
