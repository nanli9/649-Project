using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.SocialPlatforms;
using static UnityEngine.ParticleSystem;

public class generate_particles : MonoBehaviour
{
    const int binSpacing = 1;
    public int solveIteration;
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

    public GameObject hapticPluginObject;
    private HapticPlugin hapticScript;

    //bounding box
    public Vector2Int boundingBox_x;
    public Vector2Int boundingBox_y;
    public Vector2Int boundingBox_z;
    Vector3[] result;

    [SerializeField]
    ComputeShader computeShader;
    private ComputeBuffer meshPropertiesBuffer;  // Buffer for structured data
    private ComputeBuffer particlesBuffer;  // Buffer for particles
    private ComputeBuffer unsortedParticlesBuffer;  // Buffer for particles
    private ComputeBuffer globalHistogramBuffer;  // Buffer for particles
    private ComputeBuffer perBlockSumBuffer;  // Buffer for particles
    private ComputeBuffer globalHistogramPrefixSumBuffer;  // Buffer for particles
    private ComputeBuffer hapticOutputBuffer;  // Buffer for hapticOutput force
    private ComputeBuffer binBuffer;  // Buffer for histogram calculation
    private ComputeBuffer argsBuffer;
    private Bounds bounds;
    private int population;
    private float inverseRestDensity;
    private int numOfBins;
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
        public int flag;
        public Vector3 acceleration;
        float density;
        public Vector3 prePosition;
        public int padding;
        public Vector3 delta_p;
        public int padding2;
        
    };
    private void Setup()
    {
        // Boundary surrounding the meshes we will be drawing.  Used for occlusion.
        bounds = new Bounds(transform.position, Vector3.one * (range + 1));
        inverseRestDensity = 1.0f / rest_density;
        InitializeBuffers();

        hapticScript = hapticPluginObject.GetComponent<HapticPlugin>();
    }
    //initialize buffer for compute shader and indirect call
    private void InitializeBuffers()
    {
        population = init_pos_range.x * init_pos_range.y * init_pos_range.z;
        result = new Vector3[1];
        numOfBins = (boundingBox_x[1] - boundingBox_x[0]) / binSpacing *
                        (boundingBox_y[1] - boundingBox_y[0]) / binSpacing *
                        (boundingBox_z[1] - boundingBox_z[0]) / binSpacing + 1;
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
        Vector3[] hapticOutput = new Vector3[1];
        int[] binCounts = new int[numOfBins];
        Array.Clear(binCounts, 0, binCounts.Length);
        //initialize the position of the particles 
        for (int i = 0;i < init_pos_range.x; i++)
        {
            for(int j = 0; j < init_pos_range.y; j++)
            {
                for (int k = 0; k < init_pos_range.z; k++)
                {
                    particle p = new particle();
                    int index = init_pos_range.y * init_pos_range.z * i + init_pos_range.z * j + k;
                    float pos_x = 3 + i * spacing;
                    float pos_y = 0 + j * spacing;
                    float pos_z = 3 + k * spacing;
                    p.position = new Vector3(pos_x, pos_y, pos_z);
                    p.velocity = new Vector3(0,0.0f,0);
                    p.acceleration = new Vector3(0,-10.0f,0);
                    p.lambda = 0.0f;
                    particles[index] = p;


                    MeshProperties props = new MeshProperties();
                    Vector3 position = new Vector3(pos_x, pos_y, pos_z);
                    Quaternion rotation = Quaternion.identity;
                    Vector3 scale = 0.2f * Vector3.one;
                    props.mat = Matrix4x4.TRS(position, rotation, scale);
                    props.color = Color.red;
                    properties[index] = props;
                }
            }
        }
        
        meshPropertiesBuffer = new ComputeBuffer(population, MeshProperties.Size());
        particlesBuffer = new ComputeBuffer(population, 5 * 4 * sizeof(float));
        unsortedParticlesBuffer = new ComputeBuffer(population, 5 * 4 * sizeof(float));

        globalHistogramBuffer = new ComputeBuffer(numOfBins, sizeof(int));
        perBlockSumBuffer = new ComputeBuffer(numOfBins, sizeof(int));
        globalHistogramPrefixSumBuffer = new ComputeBuffer(numOfBins, sizeof(int));
        hapticOutputBuffer = new ComputeBuffer(1, 3 * sizeof(float));

        binBuffer = new ComputeBuffer(1, 3 * sizeof(float));

        meshPropertiesBuffer.SetData(properties);
        particlesBuffer.SetData(particles);
        unsortedParticlesBuffer.SetData(particles);
        hapticOutputBuffer.SetData(hapticOutput);

        globalHistogramBuffer.SetData(binCounts);
        globalHistogramPrefixSumBuffer.SetData(binCounts);
        perBlockSumBuffer.SetData(binCounts);

        instanceMaterial.SetBuffer("_Properties", meshPropertiesBuffer);
    }

    
    void Start()
    {
        Setup();
    }

    // Update is called once per frame
    void Update()
    {
        //retrieve the data from haptic device
        //hapticInteractionPoint = hapticScript.Mapped_HPos;
        //hapticInForce = hapticScript.forceInput;

        //dispatch kernels
        int Pre_solve_kernel = computeShader.FindKernel("Pre_solve");
        int Post_solve_kernel = computeShader.FindKernel("Post_solve");
        int Calculate_lambda_kernel = computeShader.FindKernel("Calculate_lambda");
        int Calculate_delta_p_kernel = computeShader.FindKernel("Calculate_delta_p");
        int Update_constrain_pos_kernel = computeShader.FindKernel("Update_constrain_pos");
        int Calculate_f_pressure_kernel = computeShader.FindKernel("Calculate_f_pressure");
        int Apply_f_pressure_kernel = computeShader.FindKernel("Apply_f_pressure");
        int Histogram_kernel = computeShader.FindKernel("Histogram");
        int Scan_kernel = computeShader.FindKernel("Scan");
        int Scan_combine_kernel = computeShader.FindKernel("Scan_per_block_combine");
        int Scatter_kernel = computeShader.FindKernel("Scatter");


        computeShader.SetFloat("deltaTime",0.016f);
        computeShader.SetInt("population", population);
        computeShader.SetFloat("inverseRestDensity", inverseRestDensity);
        //computeShader.SetVector("acceleration",new Vector4(0,-10,0,0));
        computeShader.SetVector("hapticInteractionPoint", hapticInteractionPoint);
        computeShader.SetVector("hapticIntPutForce", hapticInForce);
        computeShader.SetFloat("spacing", 1.0f);
        computeShader.SetInt("min_x", boundingBox_x.x);
        computeShader.SetInt("max_x", boundingBox_x.y);
        computeShader.SetInt("min_y", boundingBox_y.x);
        computeShader.SetInt("max_y", boundingBox_y.y);
        computeShader.SetInt("min_z", boundingBox_z.x);
        computeShader.SetInt("max_z", boundingBox_z.y);

        computeShader.SetBuffer(Pre_solve_kernel, "particles", unsortedParticlesBuffer);
        computeShader.Dispatch(Pre_solve_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        //run a histogram to do the bin couns
        computeShader.SetBuffer(Histogram_kernel, "unsortedParticles", unsortedParticlesBuffer);
        computeShader.SetBuffer(Histogram_kernel, "globalHistogram", globalHistogramBuffer);
        computeShader.Dispatch(Histogram_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        //scan to get a per workgroup size prefix sum
        computeShader.SetBuffer(Scan_kernel, "scanInput", globalHistogramBuffer);
        computeShader.SetBuffer(Scan_kernel, "scanBlockSum", perBlockSumBuffer);
        computeShader.Dispatch(Scan_kernel, Mathf.CeilToInt(numOfBins / 256f), 1, 1);

        //another scan to get the prefix sum of the last sum per work group 
        computeShader.SetBuffer(Scan_kernel, "scanInput", perBlockSumBuffer);
        computeShader.SetBuffer(Scan_kernel, "scanBlockSum", hapticOutputBuffer);
        computeShader.Dispatch(Scan_kernel, Mathf.CeilToInt(numOfBins / (256f * 256f)), 1, 1);

        //combine per work group prefix sum to global prefix sum
        computeShader.SetBuffer(Scan_combine_kernel, "globalHistogram", globalHistogramBuffer);
        computeShader.SetBuffer(Scan_combine_kernel, "scanBlockSum", perBlockSumBuffer);
        computeShader.Dispatch(Scan_combine_kernel, Mathf.CeilToInt(numOfBins / 256f), 1, 1);

        //sorted the particles scatter part
        computeShader.SetBuffer(Scatter_kernel, "particles", particlesBuffer);
        computeShader.SetBuffer(Scatter_kernel, "unsortedParticles", unsortedParticlesBuffer);
        computeShader.SetBuffer(Scatter_kernel, "globalHistogram", globalHistogramBuffer);
        //computeShader.SetBuffer(Scatter_kernel, "globalHistogramPrefixSum", globalHistogramPrefixSumBuffer);
        computeShader.Dispatch(Scatter_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        /*particle[] unsortedParticles = new particle[population];
        particle[] sortedParticles = new particle[population];
        unsortedParticlesBuffer.GetData(unsortedParticles);
        particlesBuffer.GetData(sortedParticles);

        Array.Sort(unsortedParticles, (p1, p2) => p1.flag.CompareTo(p2.flag));
        float error = 0;
        for(int i = 0; i < population; i++)
        {
            error += (unsortedParticles[i].position - sortedParticles[i].position).magnitude;
        }
        Debug.Log(error);*/
        for(int i = 0; i < solveIteration; i++)
        {
            computeShader.SetBuffer(Calculate_lambda_kernel, "globalHistogram", globalHistogramBuffer);
            computeShader.SetBuffer(Calculate_lambda_kernel, "particles", particlesBuffer);
            computeShader.Dispatch(Calculate_lambda_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

            computeShader.SetBuffer(Calculate_delta_p_kernel, "globalHistogram", globalHistogramBuffer);
            computeShader.SetBuffer(Calculate_delta_p_kernel, "particles", particlesBuffer);
            computeShader.Dispatch(Calculate_delta_p_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

            computeShader.SetBuffer(Update_constrain_pos_kernel, "particles", particlesBuffer);
            computeShader.Dispatch(Update_constrain_pos_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        }
        computeShader.SetBuffer(Post_solve_kernel, "particles", particlesBuffer);
        computeShader.SetBuffer(Post_solve_kernel, "unsortedParticles", unsortedParticlesBuffer);
        computeShader.SetBuffer(Post_solve_kernel, "transformation", meshPropertiesBuffer);
        computeShader.SetBuffer(Post_solve_kernel, "globalHistogram", globalHistogramBuffer);
        computeShader.Dispatch(Post_solve_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        
        //haptic interaction part
        computeShader.SetBuffer(Calculate_f_pressure_kernel, "particles", particlesBuffer);
        computeShader.SetBuffer(Calculate_f_pressure_kernel, "hapticOutputForce", hapticOutputBuffer);
        computeShader.Dispatch(Calculate_f_pressure_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        computeShader.SetBuffer(Apply_f_pressure_kernel, "particles", particlesBuffer);
        computeShader.Dispatch(Apply_f_pressure_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        //hapticOutputBuffer.GetData(result);
        //hapticPressureForce = result[0];
        //result[0] = new Vector3(0, 0, 0);
        //Debug.Log(result[0]);
        //find the f_pressure
        //computeShader.Dispatch(Calculate_f_pressure_kernel, Mathf.CeilToInt(population / 256f), 1, 1);


        //indirect call to render
        instanceMaterial.SetBuffer("_Properties", meshPropertiesBuffer);
        Graphics.DrawMeshInstancedIndirect(sphereMesh, 0, instanceMaterial, bounds, argsBuffer);
    }
}
