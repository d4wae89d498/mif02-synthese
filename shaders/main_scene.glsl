// Modeling - 2024.09.15
// Eric Galin
#include "_common.glsl"

struct Ray{
    vec3 o;// Origin
    vec3 d;// Direction
};

struct Material {
    vec3 color;         // Couleur du matériau
    float ambient;      // Coefficient de réflexion ambiante
    float diffuse;      // Coefficient de réflexion diffuse
    float specular;     // Coefficient de réflexion spéculaire
    float reflectivity; // Facteur de réflexion (0 = mat, 1 = miroir)
};



struct Val {
  float v; // Signed distance
  int c; // Cost
  Material mat;  // Matériau au point donné
};

// Compute point on ray
// ray : The ray
//   t : Distance
vec3 Point(Ray ray,float t)
{
    return ray.o+t*ray.d;
}


Material whiteMat = Material(vec3(1.0, 1.0, 1.0), 0.2, 0.5, 0.1, 0.0);

Material blueReflective = Material(vec3(1.0, 1.0, 1.0), 0.5, .2, 10.0, .8);

Material redMat = Material(vec3(1.0, 0.0, 0.0),  0.2, 0.5, 0.1, 0.0);


// Camera -------------------------------------------------------------------------------

// Rotation matrix around z axis
// a : Angle
mat3 Rz(float a)
{
  float sa=sin(a);float ca=cos(a);
  return mat3(ca,sa,0.,-sa,ca,0.,0.,0.,1.);
}

// Compute the ray
//      m : Mouse position
//      p : Pixel
Ray CreateRay(vec2 m,vec2 p)
{
  float a=3.*3.14*m.x;
  float le=3.5;

  // Origin
  vec3 ro=vec3(37.,0.,15.);
  ro*=Rz(a);

  // Target point
  vec3 ta=vec3(0.,0.,1.);

  // Orthonormal frame
  vec3 w=normalize(ta-ro);
  vec3 u=normalize(cross(w,vec3(0.,0.,1.)));
  vec3 v=normalize(cross(u,w));
  vec3 rd=normalize(p.x*u+p.y*v+le*w);
  return Ray(ro,rd);
}


// Object transformations ------------------------------------------------------------------

// Translate a point
// p : the point to translate
// t : translation vector
vec3 Translate(vec3 p, vec3 t) {
    return p - t;
}

// Scale (Homothety)
// p : the point to scale
// s : scaling factor (can be a vector for non-uniform scaling)
vec3 Scale(vec3 p, vec3 s) {
    return p / s;
}

// Rotate around X axis
// p : the point to rotate
// a : angle (in radians) to rotate around the X-axis
vec3 RotateX(vec3 p, float a) {
    float sa = sin(a);
    float ca = cos(a);
    return vec3(p.x, ca * p.y - sa * p.z, sa * p.y + ca * p.z);
}

// Rotate around Y axis
// p : the point to rotate
// a : angle (in radians) to rotate around the Y-axis
vec3 RotateY(vec3 p, float a) {
    float sa = sin(a);
    float ca = cos(a);
    return vec3(ca * p.x + sa * p.z, p.y, -sa * p.x + ca * p.z);
}

// Rotate around Z axis
// p : the point to rotate
// a : angle (in radians) to rotate around the Z-axis
vec3 RotateZ(vec3 p, float a) {
    float sa = sin(a);
    float ca = cos(a);
    return vec3(ca * p.x - sa * p.y, sa * p.x + ca * p.y, p.z);
}

// Transform a point with translation, scaling (homothety), and rotation
// p : the point to transform
// t : translation vector
// s : scaling vector (non-uniform scaling possible)
// r : rotation angles for X, Y, and Z axes (in radians)
vec3 Transform(vec3 p, vec3 t, vec3 s, vec3 r) {
    p = Translate(p, t);
    p = Scale(p, s);
    p = RotateX(p, r.x);
    p = RotateY(p, r.y);
    p = RotateZ(p, r.z);
    return p;
}

// Operators -------------------------------------------------------------------------------

// Union de deux primitives
Val Union(Val a, Val b) {
    // Choisir la primitive avec la distance la plus proche
    if (a.v < b.v) {
        return Val(a.v, a.c + b.c + 1, a.mat);
    } else {
        return Val(b.v, a.c + b.c + 1, b.mat);
    }
}

// Difference of two signed distance functions
Val Difference(Val a, Val b) {
    return Val(max(a.v, -b.v), a.c + b.c + 1, a.mat); // Subtract b from a
}

// Intersection of two signed distance functions
Val Intersection(Val a, Val b) {
    return Val(max(a.v, b.v), a.c + b.c + 1, a.mat); // Keep the larger distance
}


// Primitives -------------------------------------------------------------------------------


// Sphere
// p : point
// c : center of skeleton
// r : radius
Val Sphere(vec3 p,vec3 c,float r, Material mat)
{
  return Val(length(p-c)-r,1, mat);
}

// Plane
// p : point
// n : Normal of plane
// o : Point on plane
Val Plane(vec3 p, vec3 n, vec3 o, Material mat)
{
    return Val(dot((p-o),n),1, mat);
}

// Ellipsoid
// p : point to test
// c : center of the ellipsoid
// r : radii along the ellipsoid axes (as a vec3)
Val Ellipsoid(vec3 p, vec3 c, vec3 r, Material mat) {
    return Val(length((p - c) / r) - 1.0, 1, mat);
}

// Box
// p : point to test
// c : center of the box
// r : half-extents (size in each direction from the center) of the box
Val Box(vec3 p, vec3 c, vec3 r, Material mat) {
    vec3 d = abs(p - c) - r;
    return Val(length(max(d, 0.0)) + min(max(d.x, max(d.y, d.z)), 0.0), 1, mat);
}

// Cylinder
// p : point to test
// c : center of the base of the cylinder
// r : radius of the cylinder
// h : height of the cylinder
Val Cylinder(vec3 p, vec3 c, float r, float h, Material mat) {
    vec2 d = abs(vec2(length(p.xz - c.xz), p.y - c.y)) - vec2(r, h);
    return Val(min(max(d.x, d.y), 0.0) + length(max(d, 0.0)), 5, mat);
}

// Capsule
// p : point to test
// a : starting point of the capsule's line segment
// b : ending point of the capsule's line segment
// r : radius of the capsule
Val Capsule(vec3 p, vec3 a, vec3 b, float r, Material mat) {
    vec3 pa = p - a, ba = b - a;
    float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 2.0);
    return Val(length(pa - ba * h) - r, 1, mat);
}

// Torus
// p : point to test
// c : center of the torus
// t : 2D vector representing the major (t.x) and minor (t.y) radii of the torus
Val Torus(vec3 p, vec3 c, vec2 t, Material mat) {
    vec3 q = p - c;
    return Val(length(vec2(length(q.xz) - t.x, q.y)) - t.y, 3, mat);
}


// Transformed primitives ------------------------------------------------------------

// Transformed Cylinder
// p : the point in space
// c : center of the base of the cylinder
// r : radius of the cylinder
// h : height of the cylinder
// t : translation vector
// s : scaling vector
// rAngles : rotation angles (x, y, z) in radians
Val TransformedCylinder(vec3 p, vec3 c, float r, float h, vec3 t, vec3 s, vec3 rAngles) {
    p = Transform(p, t, s, rAngles);
    c = Transform(c, t, s, rAngles);
    return Cylinder(p, c, r, h, whiteMat);
}


// Primitives deformation ------------------------------------------------------------


float Noise(vec3 p) {
    return fract(sin(dot(p, vec3(12.9898, 78.233, 151.7182))) * 43758.5453);
}

float Turbulence(vec3 p, float scale) {
    float value = 0.0;
    float frequency = scale;
    for (int i = 0; i < 5; i++) {
        value += Noise(p * frequency) / frequency;
        frequency *= 2.0; // Increase frequency
    }
    return value;
}


Val DeformedPlane(vec3 p, vec3 n, vec3 o) {
    Val planeVal = Plane(p, n, o, whiteMat);
    float perturbation = Turbulence(p, 4096.); // Ajuster l'échelle du bruit
    return Val(planeVal.v + perturbation, planeVal.c, planeVal.mat);
}


Val DeformWithNoise(Val original, vec3 p) {
    // Appliquer du bruit pour déformer la distance
    float noiseScale = 0.1; // Échelle du bruit
    float noise = Noise(p * noiseScale); // Obtenir le bruit
    float deformation = 10. * noise; // Amplitude de la déformation (ajuster selon les besoins)

    return Val(original.v + deformation, original.c, original.mat);
}


// cost : 40
vec3 DeformedSphere(vec3 p, float radius) {
    float deformation = Noise(p * 2.0) * 0.5; // Ajuste l'amplitude de la déformation
    float distance = length(p) - (radius + deformation);
    return normalize(p) * (radius + deformation);
}


// Potential field of the object
// p : point
Val object(vec3 p)
{
    Val ellipsoid = Ellipsoid(p, vec3(0.0, 2.0, 2.0), vec3(2.0, 3.0, 4.0), whiteMat);
    Val box = Box(p, vec3(3.0, -8.0, -1.0), vec3(1.0, 1.0, 1.0), whiteMat);
    Val sphere = Sphere(p, vec3(-4.0, 6.0, -2.0), 4.0, blueReflective);
    Val sphere2 = Sphere(p, vec3(-9.0, 12.0, -2.0), 2.0, redMat);


    Val bsphere = Sphere(p, vec3(4.0, 6.0, -2.0), 3.0, whiteMat);

    Val cylinder = Cylinder(p, vec3(7.0, 0.0, 0.0), 1.0, 2.0, whiteMat);
    Val capsule = Capsule(p, vec3(0.0, -1.0, 5.5), vec3(0.0, 3.0, 5.5), 0.5, whiteMat);
    Val torus = Torus(p, vec3(3.0, -8.0, -2.0), vec2(1.5, 0.3), whiteMat);

    vec3 translation = vec3(0.0, 0.0, 0.0);  // No translation
    vec3 scaling = vec3(1.0, 1.0, 1.0);      // No scaling
    vec3 rotation = vec3(radians(90.0), radians(0.0), 0.0); // Rotate 90 degrees around the Y-axis
    vec3 cylinderBaseCenter = vec3(7.0, 0.0, 0.0); // Base at origin
    float cylinderRadius = 1.0; // Radius of the cylinder
    float cylinderHeight = 2.0;  // Height of the cylinder

    Val vcylinder = TransformedCylinder(p, cylinderBaseCenter, cylinderRadius, cylinderHeight, translation, scaling, rotation);

      bsphere = DeformWithNoise(bsphere, p);
    Val v = Union(bsphere, sphere);
    v = Union(v, sphere2);
    v = Union(v, Intersection(box, torus));
    v = Union(v,
        Difference(capsule, ellipsoid)
    );


    v = Union(v, cylinder);
    v = Union(v, vcylinder);


    v = Union(v, Plane(p, vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, -4.0), whiteMat));

    return v;
/*
    Val v = Union(
        Ellipsoid(p, vec3(0., 2., 2.), vec3(2., 3., 4.)),
        Box(p, vec3(3., 0., -1.), vec3(1.0, 1.0, 1.0))
    );

    v = Union(v, Torus(p, vec3(0.0, 0.0, -2.0), vec2(1.0, 0.25)));

    return v;*/


}


// Random direction in a hemisphere
// seed : Integer seed, from 0 to N
//    n : Direction of the hemisphere
vec3 Hemisphere(int seed,vec3 n)
{
    float a=fract(sin(176.19*float(seed)));// Uniform randoms
    float b=fract(sin(164.19*float(seed)));

    float u=2.*3.1415*a;// Random angle
    float v=acos(2.*b-1.);// Arccosine distribution to compensate at poles

    vec3 d=vec3(cos(u)*cos(v),sin(u)*cos(v),sin(v));// Direction
    if(dot(d,n)<0.){d=-d;}// Hemisphere

    return d;
}

// Analysis of the scalar field -----------------------------------------------------------------

const int Steps=200;// Number of steps
const float Epsilon=.01;// Marching epsilon

// Object normal
// p : point
vec3 ObjectNormal(vec3 p)
{
  const float eps=.001;
  vec3 n;
  Val val=object(p);
  float v=val.v;
  n.x=object(vec3(p.x+eps,p.y,p.z)).v-v;
  n.y=object(vec3(p.x,p.y+eps,p.z)).v-v;
  n.z=object(vec3(p.x,p.y,p.z+eps)).v-v;
  return normalize(n);
}

// Trace ray using ray marching
// ray : The ray
//   e : Maximum distance
//   h : hit
//   s : Number of steps
//   c : cost
bool SphereTrace(Ray ray,float e,out float t,out int s,out int c, out Material mat)
{
  bool h=false;

  // Start at the origin
  t=0.0;
  c=0;

  for(int i=0;i<Steps;i++)
  {
    s=i;
    vec3 p=Point(ray,t);
    Val val=object(p);
    float v=val.v;
    c+=val.c;
    // Hit object
    if(v<0.)
    {
      mat = val.mat;
      h=true;
      break;
    }
    // Move along ray
    t+=max(Epsilon,v);
    // Escape marched too far away
    if(t>e)
    {
      break;
    }
  }
  return h;
}

// Lighting -------------------------------------------------------------------------------

// Background color
// ray : Ray
vec3 background(Ray ray)
{
  return mix(vec3(.45,.55,.99),vec3(.65,.69,.99),ray.d.z*.5+.5);
}

// Shadowing
// p : Point
// n : Normal
// l : Light direction
float Shadow(vec3 p,vec3 n,vec3 l)
{
  float t;
  int s;
  int c;
  Material mat;
  bool hit=SphereTrace(Ray(p+Epsilon*n,l),100.,t,s,c, mat);
  if(!hit)
  {
    return 1.;
  }
  return 0.;
}


// Compute ambient occlusion
// p : Point
// n : Normal at point
float AmbientOcclusion(vec3 p, vec3 n)
{
    float occlusion = 0.0;
    const int samples = 16;  // Nombre d'échantillons pour l'occlusion
    float maxDistance = 5.0;  // Distance maximale de l'occlusion

    for (int i = 0; i < samples; i++)
    {
        vec3 dir = Hemisphere(i, n);  // Générer une direction aléatoire dans l'hémisphère
        float t;
        int s;
        int c;
  Material mat;

        bool hit = SphereTrace(Ray(p + Epsilon * n, dir), maxDistance, t, s, c, mat);
        if (hit)
        {
            occlusion += 1.0;
        }
    }

    occlusion = 1.0 - (occlusion / float(samples));  // Normaliser la valeur d'occlusion
    return occlusion;
}

// Modulate ambient occlusion based on direct lighting
// p : Point
// n : Normal at point
// l : Light direction
float ModulateOcclusionByLight(vec3 p, vec3 n, vec3 l)
{
    // Calculate how much the point is directly lit
    float directLightIntensity = clamp(dot(n, l), 0.0, 1.0);  // Intensité de la lumière directe

    // Calculate ambient occlusion
    float ao = AmbientOcclusion(p, n);

    // Modulate occlusion based on direct light intensity
    // If the point is fully lit (directLightIntensity ~ 1), reduce the AO
    // If the point is in shadow (directLightIntensity ~ 0), keep the full AO
    float filteredAO = mix(ao, 1.0, directLightIntensity);  // Réduire AO dans les zones éclairées
    return filteredAO;
}

// Fonction simplifiée pour tracer un rayon réfléchi
vec3 TraceReflection(vec3 p, vec3 reflectionDir, float reflectivity, float e, int depth)
{
    // Limite la profondeur de récursivité pour éviter les boucles infinies
    if (depth <= 0) return vec3(0.0); // Retourne noir si on atteint la profondeur maximale


    Ray reflectionRay = Ray(p + reflectionDir * 1.10 /** 0.010*/, reflectionDir * 1.10); // Un petit offset pour éviter l'auto-intersection
    float t = 0.0;           // Distance le long du rayon
    int s = 0;              // Nombre d'étapes
    int c = 0;              // Compteur de collisions
    Material mat;           // Matériau intersecté

    // Appel de la fonction SphereTrace pour détecter l'intersection
    bool hit = SphereTrace(reflectionRay, e, t, s, c, mat); // e est la distance d'échappement



    if (hit) {
        // Si on a touché un objet, calculez la couleur réfléchie
        vec3 reflectedColor = mat.color; // Utilisez mat pour obtenir la couleur selon vos besoins
        return reflectedColor * reflectivity; // Applique le coefficient de réflexion
    } else {
        // Sinon, retourne la couleur d'arrière-plan
        return background(reflectionRay);
    }
}



// Shading and lighting
// p : Point
// n : Normal at point
// eye : Eye direction
vec3 ShadeWithAO1(vec3 p, vec3 n, Ray eye, Material mat)
{
    // Point light
    const vec3 lp = vec3(5.0, 10.0, 25.0);

    // Light direction to point light
    vec3 l = normalize(lp - p);

    // Compute ambient occlusion
    float ao = AmbientOcclusion(p, n);  // Calculer l'occlusion ambiante

    // Ambient color with occlusion
    vec3 ambient = (0.25 + 0.25 * background(Ray(p, n))) * ao;

    // Shadow computation
    float shadow = Shadow(p, n, l);

    // Phong diffuse
    vec3 diffuse = 0.35 * clamp(dot(n, l), 0.0, 1.0) * mat.color;//vec3(1.0, 1.0, 1.0);

    // Specular
    vec3 r = reflect(eye.d, n);
    vec3 specular = 0.15 * pow(clamp(dot(r, l), 0.0, 1.0), 35.0) * vec3(1.0, 1.0, 1.0);

    // Combine ambient, diffuse, and specular components with shadow and occlusion
    vec3 c = ambient + shadow * (diffuse + specular);
    return c;
}


// Shading function with ambient occlusion
// p : Point
// n : Normal at point
// eye : Ray from the camera
vec3 ShadeWithAO2(vec3 p, vec3 n, Ray eye, Material mat)
{
    const vec3 lightPos = vec3(5.0, 10.0, 25.0);  // Position de la lumière
    vec3 l = normalize(lightPos - p);  // Direction de la lumière

    // Calculate ambient occlusion modulated by light
    float ambientOcclusion = ModulateOcclusionByLight(p, n, l);

    // Compute base lighting (diffuse + specular)
    vec3 ambient = .25 + .25 * background(Ray(p, n));
    float shadow = Shadow(p, n, l);
    vec3 diffuse = .35 * clamp(dot(n, l), 0.0, 1.0) * mat.color;
    vec3 r = reflect(eye.d, n);
    vec3 specular = .15 * pow(clamp(dot(r, l), 0.0, 1.0), 35.0) * vec3(1.0, 1.0, 1.0);

    // Combine all with the modulated occlusion
    vec3 color = ambient * ambientOcclusion + shadow * (diffuse + specular);

    return color;
}

// Fonction de shading avec Ambient Occlusion et matériaux
vec3 ShadeWithAO3(vec3 p, vec3 n, Ray eye, Material mat)
{
   // Set up the moving light position
    const float speed = .005; // Speed factor
    vec3 lightPos = vec3(50.0 * sin(iTime * speed), 10.0, 100.0 * cos(iTime * speed));

    vec3 l = normalize(lightPos - p);  // Direction de la lumière

    // Calcul de l'occlusion ambiante modulée par la lumière
    float ambientOcclusion = ModulateOcclusionByLight(p, n, l);

    // Lumière ambiante avec les coefficients de matériau
    vec3 ambient = mat.ambient * mat.color;//

    // Calcul de l'ombre (si la lumière est bloquée)
    float shadow = Shadow(p, n, l);

    // Calcul de l'éclairage diffus
    vec3 diffuse = mat.diffuse * clamp(dot(n, l), 0.0, 1.0) * mat.color;

    // Calcul de la réflexion spéculaire
    vec3 r = reflect(eye.d, n);
    vec3 specular = mat.specular * pow(clamp(dot(r, l), 0.0, 1.0), 35.0) * vec3(1.0, 1.0, 1.0);

    // Combine l'ambiant, le diffus et le spéculaire avec l'occlusion ambiante
    vec3 color = (ambient * ambientOcclusion) + shadow * (diffuse + specular);

    if (mat.reflectivity > 0.0) {
        vec3 reflectionDir = reflect(eye.d, n) * .1;
        vec3 reflection = TraceReflection(p, reflectionDir, mat.reflectivity, 1000.0, 3); // Profondeur max de 3
        color = mix(color, reflection, mat.reflectivity); // Combine la couleur de base avec la réflexion
    }


    return color;
}


// Shading and lighting
//   p : Point
//   n : Normal at point
// eye : Eye direction
vec3 Shade(vec3 p,vec3 n,Ray eye, Material mat)
{
  // Point light
  const vec3 lp=vec3(5.,10.,25.);

  // Light direction to point light
  vec3 l=normalize(lp-p);

  // Ambient color
  vec3 ambient=.25+.25*background(Ray(p,n));

  // Shadow computation
  float shadow=Shadow(p,n,l);

  // Phong diffuse
  vec3 diffuse=.35*clamp(dot(n,l),0.,1.)* mat.color;

  // Specular
  vec3 r=reflect(eye.d,n);
  vec3 specular=.15*pow(clamp(dot(r,l),0.,1.),35.)*vec3(1.,1.,1.);
  vec3 c=ambient+shadow*(diffuse+specular);
  return c;
}

// Shading according to the number of steps in sphere tracing
// n : Number of steps
vec3 ShadeSteps(int n,int m)
{
  float t=float(n)/(float(m));
  return.5+mix(vec3(.05,.05,.5),vec3(.65,.39,.65),t);
}


// Image
void mainImage(out vec4 color,in vec2 pxy)
{
  // Convert pixel coordinates
    vec2 pixel=(-iResolution.xy+2.*pxy)/iResolution.y;

  // Mouse
  vec2 m=iMouse.xy/iResolution.xy;

  // Camera
  Ray ray=CreateRay(m,pixel);

  // Trace ray

  // Hit and number of steps
  float t=0.0;
  int s=0;

  int c;
  Material mat;
  bool hit=SphereTrace(ray,100.,t,s, c, mat);

  // Shade background
  vec3 rgb=background(ray);

  if(hit)
  {
    // Position
    vec3 p=Point(ray,t);

    // Compute normal
    vec3 n=ObjectNormal(p);

    // Shade object with light
    rgb=ShadeWithAO3(p,n,ray, mat);


  }

  // Uncomment this line to shade image with false colors representing the number of steps
  //  rgb=ShadeSteps(s,Steps);

   // Uncomment this line to shade cost
  //rgb=ShadeSteps(c,500);
  color=vec4(rgb,1.);
}
