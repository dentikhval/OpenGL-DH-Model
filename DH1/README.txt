This is a 3D visualisation of a robot defined by Denavit-Hartenberg (DH) parameters. This program draws it on the screen, and you can control individual joints by pressing keys on your keyboard. XYZ end effector position is conveniently displayed in the console output.
It can be a good way to verify if your input data for a robot actually works as you expected and the robot defined by DH params actually looks like your real life robot.
For an elevated experience, you can substitute the STL models attached with the model of your robot - just for fun.

The author decided by some reason to hardcode the robot dimensions in function void Robot::configurarTH(), so you will have to go hunting for them manually and guess where to set a, d, alpha and theta. Why not simply implement a 4 x 6 matrix as input?!
Console output seems to be XYZ end effector position - nice!

I know there are better programs now like ROS, Moveit, openRave etc., but they are too long to learn and have their own special quirks. This, though, is simple, fully configurable and runs locally.

Why do all experianced developers think everyone will automatically understand what to do with their code?!!!
I had no idea about it, but some Arduino programming skills.

First, I had no idea which libraries to include as compile time inputs, because the original authot did not bother to provide them.
Second, there was some issue with file linking, and the main file couldn't find the other files.
I can't care less what was causing it, I just want a working simulation!
So I just copied all files in the correct order into a single one, and now it works fine! Hehe...
Also, the author used a numpad for control... Who even has a numpad today?! I was running this from my laptop, had to remap skeys in a more logical arrangement.

To run this simulation, if you know nothing about OpenGL and have never run a single C++ program in Windows:

(1) Follow the attached PDF guide completely until you can draw a test triangle
(1.1) Do not forget to add the folder C:\mingw64\x86_64-w64-mingw32 to your system PATH environment variable. Not USER, SYSTEM!
(1.2) You will need to install both glut and glew libraries, don't skip that

(2) Open Command prompt or Powershell AS ADMIN and go to the directory where this file and main.cpp is:
cd C:\mingw64\x86_64-w64-mingw32\DH1

(3) Run the following command:

 g++ main.cpp -o DH1  -lopengl32  -lglu32 -lgdi32

It will create an EXE file which you can run and control the robot with your keyboard.
The compiler will give you some errors like this:
 warning: control reaches end of non-void function [-Wreturn-type]
Don't worry, these do not prevent you from successful compiling.

Joint number: key for positive direction, key for negative direction. Do not press "+" or "-"
j1: r+ f-
j2: t+ g-
j3: y+ h-
j4: u+ j-
j5: i+ k-
j6: o+ l-
j7(gripper) rotation: p- (negative only)

Camera movement: W,A,S,D

Exit: esc key