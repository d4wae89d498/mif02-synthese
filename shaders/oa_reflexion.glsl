// Modeling - 2024.09.15
// Eric Galin
#include "_common.glsl"
#line 5

struct Ray {
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

struct Light {
    vec3 position;      // Position for point lights
    vec3 direction;     // Direction for directional lights
    vec3 color;         // Light color
    float intensity;    // Light intensity
    bool isDirectional; // Whether the light is directional
};

const float sunSpeed = .10; // Speed factor
vec3 sunPos = vec3(50.0 * sin(iTime * sunSpeed), 10.0, 100.0 * cos(iTime * sunSpeed));

const int numLights = 3; // Nombre de lumières
Light lights[numLights];


void initLights() {

    lights[0] = Light(vec3(10.0, 10.0, 10.0), vec3(0.0), vec3(1.0, 0.9, 0.8), 1.0, false);  // Lumière ponctuelle
    lights[1] = Light(vec3(0.0), normalize(vec3(-1.0, -1.0, -0.5)), vec3(0.8, 0.8, 1.0), 0.7, true);  // Lumière directionnelle
    lights[2] = Light(vec3(-10.0, 5.0, -10.0), vec3(0.0), vec3(0.7, 1.0, 0.7), 0.8, false);  // Deuxième lumière ponctuelle
}



// Compute point on ray
// ray : The ray
//   t : Distance
vec3 RayPoint(Ray ray,float t)
{
    return ray.o+t*ray.d;
}

const float defaultAmbiantFactor = .8;

Material concreteMat = Material(vec3(0.6, 0.6, 0.6), defaultAmbiantFactor, 0.5, 0.1, 0.0);

Material greenMat = Material(vec3(0.1882, 1.0, 0.4471), defaultAmbiantFactor, 0.7, 0.3, 0.0);

Material cyanMat = Material(vec3(0.0, 0.898, 1.0), defaultAmbiantFactor, 0.7, 0.3, 0.0);

Material orangeMat = Material(vec3(1.0, 0.749, 0.0), defaultAmbiantFactor, 0.7, 0.3, 0.0);

Material purpleMat = Material(vec3(1.0, 0.0, 0.8), defaultAmbiantFactor, 0.7, 0.3, 0.0);

Material mirrorMat = Material(vec3(1.0, 1.0, 1.0), defaultAmbiantFactor, .2, 10.0, .9);

Material rubyMat = Material(vec3(0.9, 0.1, 0.1), defaultAmbiantFactor, 0.7, 0.9, 0.2);

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

// Main object ---------------------------------------------------------------------

// Potential field of the object
// p : point
Val object(vec3 p)
{
    Val ellipsoid = Ellipsoid(p, vec3(0.0, 2.0, 2.0), vec3(2.0, 3.0, 4.0), concreteMat);
    Val box = Box(p, vec3(3.0, -8.0, -1.0), vec3(1.0, 1.0, 1.0), purpleMat);
    Val sphere = Sphere(p, vec3(4.0, -4.0, -2.0), 3.0, mirrorMat);
    Val sphere2 = Sphere(p, vec3(-2.0, 3.0, -2.0), 2.0, rubyMat);

    Val bsphere = Sphere(p, vec3(4.0, 6.0, -2.0), 1.0, greenMat);

    Val cylinder = Cylinder(p, vec3(7.0, 0.0, 0.0), 1.0, 2.0, orangeMat);
    Val capsule = Capsule(p, vec3(0.0, -1.0, 5.5), vec3(0.0, 3.0, 5.5), 0.5, cyanMat);
    Val torus = Torus(p, vec3(3.0, -8.0, -2.0), vec2(1.5, 0.3), purpleMat);

	Val vcylinder = Cylinder(RotateX(p, radians(90.)), vec3(7.0, 0.0, 0.0), 1.0, 2.0, orangeMat);

    Val v = Union(bsphere, sphere);
    v = Union(v, sphere2);
    v = Union(v, Intersection(box, torus));
    v = Union(v,
        Difference(capsule, ellipsoid)
    );


    v = Union(v, cylinder);
    v = Union(v, vcylinder);


    v = Union(v, Plane(p, vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, -4.0), concreteMat));

	Val mirror = Box(
		RotateX(RotateZ(p, radians(90.)), radians(-4.)),
		vec3(-2., -7.0, -3.0), vec3(13.0, 0.1, 12.0), mirrorMat);
	//mirror.RotateY(radians(10), mirror);
    v = Union(
		v,
		mirror//Box(p, vec3(3.0, -13.0, -4.0), vec3(10.0, 0.1, 10.0), mirrorMat)//, rotationMatrix(radians(-6.0), vec3(1.0, 0.0, 0.0)))
		);

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
    vec3 p=RayPoint(ray,t);
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

vec3 Shade(vec3 p, vec3 n, Ray eye, Material mat)
{
    // Light direction to point light
    vec3 l = normalize(sunPos - p);

    // Ambient component based on material properties
    vec3 ambient = mat.ambient * mat.color;  // Material ambient color scaled by ambient strength

    // Shadow computation
    float shadow = Shadow(p, n, l);

    // Phong diffuse component based on material's diffuse property
    vec3 diffuse = mat.diffuse * clamp(dot(n, l), 0.0, 1.0) * mat.color;

    // Specular component based on material's specular property
    vec3 r = reflect(eye.d, n);
    vec3 specular = mat.specular * pow(clamp(dot(r, l), 0.0, 1.0), 35.0) * vec3(1.0);

    // Combine all components: ambient, diffuse, specular, shadow
    vec3 color = ambient + shadow * (diffuse + specular);

    return color;
}


// Compute ambient occlusion
// p : Point
// n : Normal at point
float AmbientOcclusion(vec3 p, vec3 n)
{
    float occlusion = 0.0;
    const int samples = 128;  // Nombre d'échantillons pour l'occlusion
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

// Shading and lighting
// p : Point
// n : Normal at point
// eye : Eye direction
vec3 ShadeWithAO1(vec3 p, vec3 n, Ray eye, Material mat)
{
    // Light direction to point light (using the same sunPos from your original code)
    vec3 l = normalize(sunPos - p);

    // Compute ambient occlusion
    float ao = AmbientOcclusion(p, n);  // Compute the ambient occlusion

    // Call the modified Shade function (now using material properties)
    vec3 color = Shade(p, n, eye, mat);  // This will return color with the lighting components

    // Combine the result from Shade with ambient occlusion and shadow
    color = color * ao;  // Apply ambient occlusion by scaling the color

    return color;
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
    float filteredAO = mix(ao, 1., directLightIntensity);  // Réduire AO dans les zones éclairées
    return filteredAO;
}

// Helper function to compute lighting based on material and lighting
vec3 ShadeWithAO2(vec3 p, vec3 n, Ray eye, Material mat)
{
	vec3 l = normalize(sunPos - p);
    float ao = ModulateOcclusionByLight(p, n, l);

     // Call the modified Shade function (now using material properties)
    vec3 color = Shade(p, n, eye, mat);  // This will return color with the lighting components

    // Combine the result from Shade with ambient occlusion and shadow
    color = color * ao;  // Apply ambient occlusion by scaling the color

    return color;
}

vec3 TraceReflection(vec3 p, vec3 reflectionDir, float reflectivity, float e)
{
    Ray reflectionRay = Ray(p + reflectionDir, reflectionDir);
    float t = 0.0;           // Distance le long du rayon
    int s = 0;              // Nombre d'étapes
    int c = 0;              // Compteur de collisions
    Material mat;           // Matériau intersecté

    bool hit = SphereTrace(reflectionRay, e, t, s, c, mat); // e est la distance d'échappement
    if (hit) {
        // Si on a touché un objet, calculez la couleur réfléchie
        vec3 reflectedColor = mat.color;

        return reflectedColor * reflectivity;
    } else {
        // Sinon, retourne la couleur d'arrière-plan
        return background(reflectionRay);
    }
}

vec3 TraceReflectionWithShadows(vec3 p, vec3 reflectionDir, float reflectivity, float e)
{
    Ray reflectionRay = Ray(p + reflectionDir, reflectionDir);
    float t = 0.0;           // Distance le long du rayon
    int s = 0;              // Nombre d'étapes
    int c = 0;              // Compteur de collisions
    Material mat;           // Matériau intersecté

    bool hit = SphereTrace(reflectionRay, e, t, s, c, mat); // e est la distance d'échappement
    if (hit) {
        // Si on a touché un objet, calculez la couleur réfléchie
        vec3 reflectedColor = mat.color;

		vec3 hitPos = RayPoint(reflectionRay, t);
        vec3 hitNormal = ObjectNormal(hitPos);
        reflectedColor = ShadeWithAO2(hitPos, hitNormal, reflectionRay, mat);

		return reflectedColor * reflectivity;
       // return reflectedColor * reflectivity;
    } else {
        // Sinon, retourne la couleur d'arrière-plan
        return background(reflectionRay);
    }
}


// Fonction de shading avec Ambient Occlusion et matériaux
vec3 ShadeWithAO2AndReflection(vec3 p, vec3 n, Ray eye, Material mat)
{


    // Combine l'ambiant, le diffus et le spéculaire avec l'occlusion ambiante
    vec3 color = ShadeWithAO2(p, n, eye, mat);

    if (mat.reflectivity > 0.0) {
        vec3 reflectionDir = normalize(reflect(eye.d, n) );
        vec3 reflection = TraceReflection(p, reflectionDir, mat.reflectivity, 1000.0);
        color = mix(color, reflection, mat.reflectivity); // Combine la couleur de base avec la réflexion
    }


    return color;
}

// Fonction de shading avec Ambient Occlusion et matériaux
vec3 ShadeWithAO2AndReflectionWithShadows(vec3 p, vec3 n, Ray eye, Material mat)
{


    // Combine l'ambiant, le diffus et le spéculaire avec l'occlusion ambiante
    vec3 color = ShadeWithAO2(p, n, eye, mat);

    if (mat.reflectivity > 0.0) {
        vec3 reflectionDir = normalize(reflect(eye.d, n) );
        vec3 reflection = TraceReflectionWithShadows(p, reflectionDir, mat.reflectivity, 1000.0);
        color = mix(color, reflection, mat.reflectivity); // Combine la couleur de base avec la réflexion
    }


    return color;
}

const int maxReflections = 5;


vec3 ShadeWithAO2AndNestedReflection(vec3 p, vec3 n, Ray eye, Material mat)
{

    vec3 color = ShadeWithAO2(p, n, eye, mat);

    if (mat.reflectivity > 0.0) {
        vec3 reflectionDir = normalize(reflect(eye.d, n));
        float reflectivity = mat.reflectivity;

        for (int i = 0; i < maxReflections; i++) {
			// if (reflectivity <= 0.0001)
			//	break;
            Ray reflectionRay = Ray(p + reflectionDir, reflectionDir);
            float t = 0.0;
            int s = 0, c = 0;
            Material matHit;

            bool hit = SphereTrace(reflectionRay, 1000.0, t, s, c, matHit);
              vec3 hitPos = RayPoint(reflectionRay, t);
                vec3 hitNormal = ObjectNormal(hitPos);

			if (hit) {

    			vec3 reflectedColor = ShadeWithAO2(hitPos, hitNormal, reflectionRay, matHit);


                // Mix with accumulated color
                color = mix(color, reflectedColor, reflectivity);

                // Prepare for the next iteration
                reflectivity *= matHit.reflectivity; // Reduce reflectivity for subsequent bounces
                reflectionDir = normalize(reflect(reflectionRay.d, hitNormal)); // Update direction for next bounce
                p = hitPos; // Move position to the hit point for the next reflection bounce
            } else {
				color = mix(color, background(reflectionRay), reflectivity);

                break; // Exit if no hit (background reached)
            }
        }
    }

    return color;
}

// Image
void mainImage(out vec4 color,in vec2 pxy)
{

	initLights();

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
    vec3 p=RayPoint(ray,t);

    // Compute normal
    vec3 n=ObjectNormal(p);

    // Shade object with light
    rgb=ShadeWithAO2(p,n,ray, mat);


  }
  color=vec4(rgb,1.);
}
