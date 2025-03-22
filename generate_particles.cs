using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.SocialPlatforms;

public class generate_particles : MonoBehaviour
{
    public float range;
    public Vector3Int init_pos_range;
    public float rest_density;
    public float spacing;
    public Material instanceMaterial;
    public Mesh sphereMesh;
    //haptic
    public Vector3 hapticInteractionPoint; //read only
    public Vector3 hapticPressureForce; // what I will update
    public Vector3 hapticInForce; //read only

    //bounding box
    public Vector2Int boundingBox_x;
    public Vector2Int boundingBox_y;
    public Vector2Int boundingBox_z;

    [SerializeField]
    ComputeShader computeShader;
    private ComputeBuffer meshPropertiesBuffer;  // Buffer for structured data
    private ComputeBuffer particlesBuffer;  // Buffer for particles
    private ComputeBuffer argsBuffer;
    private Bounds bounds;
    private int population;
    private float inverseRestDensity;
    // Start is called before the first frame update
    private struct MeshProperties
    {
        public Matrix4x4 mat;
        public Vector4 color;

        public static int Size()
        {
            return
                sizeof(float) * 4 * 4 + // matrix;
                sizeof(float) * 4;      // color;
        }
    }
    private struct particle
    {
        public Vector3 velocity;
        public float lambda;
        public Vector3 position;
        float padding;
    };
    private void Setup()
    {
        // Boundary surrounding the meshes we will be drawing.  Used for occlusion.
        bounds = new Bounds(transform.position, Vector3.one * (range + 1));
        inverseRestDensity = 1.0f / rest_density;
        InitializeBuffers();
    }
    //initialize buffer for compute shader and indirect call
    private void InitializeBuffers()
    {
        population = init_pos_range.x * init_pos_range.y * init_pos_range.z;

        // Argument buffer used by DrawMeshInstancedIndirect.
        uint[] args = new uint[5] { 0, 0, 0, 0, 0 };
        // Arguments for drawing mesh.
        // 0 == number of triangle indices, 1 == population, others are only relevant if drawing submeshes.
        args[0] = (uint)sphereMesh.GetIndexCount(0);
        args[1] = (uint)population;
        args[2] = (uint)sphereMesh.GetIndexStart(0);
        args[3] = (uint)sphereMesh.GetBaseVertex(0);
        argsBuffer = new ComputeBuffer(1, args.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        argsBuffer.SetData(args);

        // Initialize buffer with the given population.
        MeshProperties[] properties = new MeshProperties[population];
        particle[] particles = new particle[population];
        //initialize the position of the particles 
        for (int i = 0;i < init_pos_range.x; i++)
        {
            for(int j = 0; j < init_pos_range.y; j++)
            {
                for (int k = 0; k < init_pos_range.z; k++)
                {
                    particle p = new particle();
                    int index = init_pos_range.y * init_pos_range.z * i + init_pos_range.z * j + k;
                    float pos_x = i * spacing;
                    float pos_y = j * spacing;
                    float pos_z = k * spacing;
                    p.position = new Vector3(pos_x, pos_y, pos_z);
                    p.velocity = new Vector3(0,0.01f,0);
                    p.lambda = 0.0f;
                    particles[index] = p;

                    MeshProperties props = new MeshProperties();
                    Vector3 position = new Vector3(pos_x, pos_y, pos_z);
                    Quaternion rotation = Quaternion.identity;
                    Vector3 scale = 0.5f * Vector3.one;
                    props.mat = Matrix4x4.TRS(position, rotation, scale);
                    props.color = Color.Lerp(Color.red, Color.blue, Random.value);
                    properties[index] = props;
                }
            }
        }

        meshPropertiesBuffer = new ComputeBuffer(population, MeshProperties.Size());
        particlesBuffer = new ComputeBuffer(population, 2 * 4 * sizeof(float));

        meshPropertiesBuffer.SetData(properties);
        particlesBuffer.SetData(particles);

        instanceMaterial.SetBuffer("_Properties", meshPropertiesBuffer);
    }

    
    void Start()
    {
        Setup();
    }

    // Update is called once per frame
    void Update()
    {
        //dispatch kernels
        int Pre_solve_kernel = computeShader.FindKernel("Pre_solve");
        int Post_solve_kernel = computeShader.FindKernel("Post_solve");
        int Calculate_lambda_kernel = computeShader.FindKernel("Calculate_lambda");
        int Calculate_delta_p_kernel = computeShader.FindKernel("Calculate_delta_p");
        computeShader.SetFloat("deltaTime",0.01f);
        computeShader.SetInt("population", population);
        computeShader.SetFloat("inverseRestDensity", inverseRestDensity);
        computeShader.SetVector("acceleration",new Vector4(0,-10,0,0));

        computeShader.SetBuffer(Pre_solve_kernel, "particles", particlesBuffer);
        computeShader.Dispatch(Pre_solve_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        //find a neighbor

        computeShader.SetBuffer(Calculate_lambda_kernel, "particles", particlesBuffer);
        computeShader.Dispatch(Calculate_lambda_kernel, Mathf.CeilToInt(population / 256f), 1, 1);


        computeShader.SetBuffer(Calculate_delta_p_kernel, "particles", particlesBuffer);
        computeShader.Dispatch(Calculate_delta_p_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        computeShader.SetBuffer(Post_solve_kernel, "particles", particlesBuffer);
        computeShader.SetBuffer(Post_solve_kernel, "transformation", meshPropertiesBuffer);
        computeShader.Dispatch(Post_solve_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        //find the f_pressure



        //indirect call to render
        instanceMaterial.SetBuffer("_Properties", meshPropertiesBuffer);
        Graphics.DrawMeshInstancedIndirect(sphereMesh, 0, instanceMaterial, bounds, argsBuffer);
        //int kernelHandle = computeShader.FindKernel("Generate_pos");

        // Set the buffer to the compute shader
        /*computeShader.SetBuffer(kernelHandle, "positions", dataBuffer);
        computeShader.Dispatch(kernelHandle, instances_x / numThread_x, instances_y / numThread_y, instances_z / numThread_z);
        instanceMaterial.SetBuffer("_InstanceDataBuffer", dataBuffer);*/

        /*Vector3[] results = new Vector3[instances_x * instances_y * instances_z];
        dataBuffer.GetData(results);*/

        // Log the results
        /*for (int i = 0; i < results.Length; i++)
        {
            Debug.Log($"Result {i}: {results[i].ToString()}");
        }
        // Release the buffer after use
        dataBuffer.Release();*/
    }
    void OnRenderObject()
    {
        instanceMaterial.SetPass(0);

        // Draw the instances
        //Graphics.DrawMeshInstancedProcedural(sphereMesh, 0, instanceMaterial, new Bounds(Vector3.zero, new Vector3(100, 100, 100)), instances_x * instances_y * instances_z);

    }
}
