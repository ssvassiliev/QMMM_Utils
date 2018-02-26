typedef struct vector {
  float x; float y; float z;
} vec3d;

vec3d veccross(vec3d, vec3d);
vec3d vecsub(vec3d, vec3d);
vec3d vecadd(vec3d, vec3d);
vec3d vecscale(vec3d, float);
float vecdot(vec3d, vec3d);
float dist3d(vec3d, vec3d);

vec3d veccross(vec3d b, vec3d c)
{
  vec3d a;
  a.x = b.y*c.z - b.z*c.y;
  a.y = b.z*c.x - b.x*c.z;
  a.z = b.x*c.y - b.y*c.x;
  return a;
}

vec3d vecsub(vec3d b, vec3d c)
{
  b.x -= c.x;
  b.y -= c.y;
  b.z -= c.z;
  return b;
}

vec3d vecadd(vec3d b, vec3d c)
{
  b.x += c.x;
  b.y += c.y;
  b.z += c.z;
  return b;
}

float vecdot(vec3d b, vec3d c)
{
  return (b.x*c.x + b.y*c.y + b.z*c.z);
}

vec3d  vecscale(vec3d b, float s)
{
 b.x*=s;
 b.y*=s;
 b.z*=s;
 return b;
}

float dist3d(vec3d b, vec3d c)
{
  return (float)sqrt((b.x-c.x)*(b.x-c.x)+(b.y-c.y)*(b.y-c.y)+(b.z-c.z)*(b.z-c.z));
}
