using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.SocialPlatforms;

public class PBF : MonoBehaviour
{
    const int binSpacing = 1;
    [Header("Simulation Parameters")]
    public int solveIteration;
    public Vector3Int init_pos_range;
    public float rest_density;
    public float init_spacing;
    public Material instanceMaterial;
    public Mesh sphereMesh;
    //haptic
    [HideInInspector] public Vector3 hapticInteractionPoint; //read only
    [HideInInspector] public Vector3 hapticPressureForce; // what I will update
    [HideInInspector] public Vector3 hapticInForce; //read only
    [HideInInspector] public Vector3 hapticInVelocity; //read only

    public GameObject hapticPluginObject;
    private HapticPlugin hapticScript;

    //bounding box
    public Vector2Int boundingBox_x = new Vector2Int(0, 30);
    public Vector2Int boundingBox_y = new Vector2Int(0, 20);
    public Vector2Int boundingBox_z = new Vector2Int(0, 20);
    float[] result;

    public ComputeShader Simulation;
    public ComputeShader Utility;

    private ComputeBuffer meshPropertiesBuffer;  // Buffer for structured data
    private ComputeBuffer particlesBuffer;  // Buffer for particles
    private ComputeBuffer unsortedParticlesBuffer;  // Buffer for particles
    private ComputeBuffer globalHistogramBuffer;  // Buffer for particles
    private ComputeBuffer perBlockSumBuffer;  // Buffer for particles
    //private ComputeBuffer globalHistogramPrefixSumBuffer;  // Buffer for particles
    private ComputeBuffer hapticOutputBuffer;  // Buffer for hapticOutput force
    private ComputeBuffer binBuffer;  // Buffer for histogram calculation
    private ComputeBuffer argsBuffer;
    private ComputeBuffer neighborListBuffer;
    private Bounds bounds;
    private int population;
    private float inverseRestDensity;
    private int numOfBins;
    private float hapticInteractionPointDensity;

    //private const int MAX_NUM_NEIGHBOR = 70;
    Vector3 densityTextureDim = new Vector3(64, 64, 64);
    [HideInInspector] public RenderTexture DensityMap;
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
        const float range = 1.0f;

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
        result = new float[3];
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
                    float pos_x = 2 + i * init_spacing;
                    float pos_y = 2 + j * init_spacing;
                    float pos_z = 0 + k * init_spacing;
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

        //neighborListBuffer = new ComputeBuffer(population * MAX_NUM_NEIGHBOR, sizeof(int));

        globalHistogramBuffer = new ComputeBuffer(numOfBins, sizeof(int));
        perBlockSumBuffer = new ComputeBuffer(numOfBins, sizeof(int));
        //globalHistogramPrefixSumBuffer = new ComputeBuffer(numOfBins, sizeof(int));
        hapticOutputBuffer = new ComputeBuffer(population, 3 * sizeof(float));

        binBuffer = new ComputeBuffer(1, 3 * sizeof(float));

        meshPropertiesBuffer.SetData(properties);
        particlesBuffer.SetData(particles);
        unsortedParticlesBuffer.SetData(particles);

        globalHistogramBuffer.SetData(binCounts);
        //globalHistogramPrefixSumBuffer.SetData(binCounts);
        perBlockSumBuffer.SetData(binCounts);

        instanceMaterial.SetBuffer("_Properties", meshPropertiesBuffer);

        RenderTextureFormat format = RenderTextureFormat.RFloat;
        DensityMap = new RenderTexture((int)densityTextureDim[0], (int)densityTextureDim[1], 0, format);
        DensityMap.dimension = UnityEngine.Rendering.TextureDimension.Tex3D;
        DensityMap.volumeDepth = (int)densityTextureDim[2];
        DensityMap.enableRandomWrite = true;
        DensityMap.useMipMap = true;

        DensityMap.Create();
    }


    void Start()
    {
        Setup();
    }

    // Update is called once per frame
    void Update()
    {
        //retrieve the data from haptic device
        hapticInteractionPoint = hapticScript.Mapped_HPos;
        hapticInForce = hapticScript.forceInput;
        hapticInVelocity = hapticScript.CurrentVelocity;
        Debug.Log(hapticInVelocity);

        //dispatch kernels
        int Pre_solve_kernel = Simulation.FindKernel("Pre_solve");
        int Post_solve_kernel = Simulation.FindKernel("Post_solve");
        int Calculate_lambda_kernel = Simulation.FindKernel("Calculate_lambda");
        int Calculate_delta_p_kernel = Simulation.FindKernel("Calculate_delta_p");
        int Update_constrain_pos_kernel = Simulation.FindKernel("Update_constrain_pos");
        int Calculate_haptic_density_kernel = Simulation.FindKernel("Calculate_haptic_density");
        int Calculate_f_pressure_kernel = Simulation.FindKernel("Calculate_f_pressure");
        int Calculate_f_viscosity_kernel = Simulation.FindKernel("Calculate_f_viscosity");
        int Apply_f_pressure_kernel = Simulation.FindKernel("Apply_f_pressure");
        int Histogram_kernel = Simulation.FindKernel("Histogram");
        int Scan_kernel = Utility.FindKernel("Scan");
        int Scan_combine_kernel = Utility.FindKernel("Scan_per_block_combine");
        int Reduce_kernel = Utility.FindKernel("Reduce");

        int Scatter_kernel = Simulation.FindKernel("Scatter");
        //int FindNeighborList_kernel = Simulation.FindKernel("FindNeighborList");

        int UpdateDensityTexture_kernel = Simulation.FindKernel("UpdateDensityTexture");


        Simulation.SetFloat("deltaTime",0.016f);
        Simulation.SetInt("population", population);
        Simulation.SetFloat("inverseRestDensity", inverseRestDensity);
        //Simulation.SetVector("acceleration",new Vector4(0,-10,0,0));
        Simulation.SetVector("hapticInteractionPoint", hapticInteractionPoint);
        Simulation.SetVector("hapticIntPutForce", hapticInForce);
        Simulation.SetVector("densityMapSize", densityTextureDim);
        Simulation.SetFloat("spacing", 1.0f);
        Simulation.SetInt("numOfParticles", population);

        //Simulation.SetInt("MAX_NUM_NEIGHBOR", MAX_NUM_NEIGHBOR);

        Simulation.SetInt("min_x", boundingBox_x[0]);
        Simulation.SetInt("max_x", boundingBox_x[1]);
        Simulation.SetInt("min_y", boundingBox_y[0]);
        Simulation.SetInt("max_y", boundingBox_y[1]);
        Simulation.SetInt("min_z", boundingBox_z[0]);
        Simulation.SetInt("max_z", boundingBox_z[1]);

        Simulation.SetInt("num_bin_x", (boundingBox_x[1] - boundingBox_x[0]) / binSpacing);
        Simulation.SetInt("num_bin_y", (boundingBox_y[1] - boundingBox_y[0]) / binSpacing);
        Simulation.SetInt("num_bin_z", (boundingBox_z[1] - boundingBox_z[0]) / binSpacing);



        Simulation.SetBuffer(Pre_solve_kernel, "particles", unsortedParticlesBuffer);
        Simulation.Dispatch(Pre_solve_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        //run a histogram to do the bin couns
        Simulation.SetBuffer(Histogram_kernel, "unsortedParticles", unsortedParticlesBuffer);
        Simulation.SetBuffer(Histogram_kernel, "globalHistogram", globalHistogramBuffer);
        Simulation.Dispatch(Histogram_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        //scan to get a per workgroup size prefix sum
        Utility.SetBuffer(Scan_kernel, "scanInput", globalHistogramBuffer);
        Utility.SetBuffer(Scan_kernel, "scanBlockSum", perBlockSumBuffer);
        Utility.Dispatch(Scan_kernel, Mathf.CeilToInt(numOfBins / 256f), 1, 1);

        //another scan to get the prefix sum of the last sum per work group 
        Utility.SetBuffer(Scan_kernel, "scanInput", perBlockSumBuffer);
        Utility.SetBuffer(Scan_kernel, "scanBlockSum", hapticOutputBuffer);
        Utility.Dispatch(Scan_kernel, Mathf.CeilToInt(numOfBins / (256f * 256f)), 1, 1);

        //combine per work group prefix sum to global prefix sum
        Utility.SetBuffer(Scan_combine_kernel, "scanPerBlockResult", globalHistogramBuffer);
        Utility.SetBuffer(Scan_combine_kernel, "scanBlockSum", perBlockSumBuffer);
        Utility.Dispatch(Scan_combine_kernel, Mathf.CeilToInt(numOfBins / 256f), 1, 1);

        //sorted the particles scatter part
        Simulation.SetBuffer(Scatter_kernel, "particles", particlesBuffer);
        Simulation.SetBuffer(Scatter_kernel, "unsortedParticles", unsortedParticlesBuffer);
        Simulation.SetBuffer(Scatter_kernel, "globalHistogram", globalHistogramBuffer);
        //Simulation.SetBuffer(Scatter_kernel, "globalHistogramPrefixSum", globalHistogramPrefixSumBuffer);
        Simulation.Dispatch(Scatter_kernel, Mathf.CeilToInt(population / 256f), 1, 1);


        //Simulation.SetBuffer(FindNeighborList_kernel, "particles", particlesBuffer);
        //Simulation.SetBuffer(FindNeighborList_kernel, "globalHistogram", globalHistogramBuffer);
        //Simulation.SetBuffer(FindNeighborList_kernel, "neighborList", neighborListBuffer);
        //Simulation.Dispatch(FindNeighborList_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

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
        for (int i = 0; i < solveIteration; i++)
        {
            Simulation.SetBuffer(Calculate_lambda_kernel, "globalHistogram", globalHistogramBuffer);
            Simulation.SetBuffer(Calculate_lambda_kernel, "particles", particlesBuffer);
            //Simulation.SetBuffer(Calculate_lambda_kernel, "neighborList", neighborListBuffer);
            Simulation.Dispatch(Calculate_lambda_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

            Simulation.SetBuffer(Calculate_delta_p_kernel, "globalHistogram", globalHistogramBuffer);
            Simulation.SetBuffer(Calculate_delta_p_kernel, "particles", particlesBuffer);
            //Simulation.SetBuffer(Calculate_delta_p_kernel, "neighborList", neighborListBuffer);
            Simulation.Dispatch(Calculate_delta_p_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

            Simulation.SetBuffer(Update_constrain_pos_kernel, "particles", particlesBuffer);
            Simulation.Dispatch(Update_constrain_pos_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        }
        Simulation.SetBuffer(Post_solve_kernel, "particles", particlesBuffer);
        Simulation.SetBuffer(Post_solve_kernel, "unsortedParticles", unsortedParticlesBuffer);
        Simulation.SetBuffer(Post_solve_kernel, "transformation", meshPropertiesBuffer);
        Simulation.Dispatch(Post_solve_kernel, Mathf.CeilToInt(population / 256f), 1, 1);


        Simulation.SetBuffer(UpdateDensityTexture_kernel, "particles", particlesBuffer);
        Simulation.SetBuffer(UpdateDensityTexture_kernel, "globalHistogram", globalHistogramBuffer);
        Simulation.SetTexture(UpdateDensityTexture_kernel, "DensityMap", DensityMap);
        Simulation.Dispatch(UpdateDensityTexture_kernel, Mathf.CeilToInt(densityTextureDim[0] / 8), Mathf.CeilToInt(densityTextureDim[1] / 8), Mathf.CeilToInt(densityTextureDim[2] / 8));

        //haptic interaction part
        Simulation.SetBuffer(Calculate_haptic_density_kernel, "particles", unsortedParticlesBuffer);
        Simulation.SetBuffer(Calculate_haptic_density_kernel, "hapticOutputForce", hapticOutputBuffer);
        Simulation.Dispatch(Calculate_haptic_density_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        //per block result
        Utility.SetBuffer(Reduce_kernel, "reduceInput", hapticOutputBuffer);
        Utility.Dispatch(Reduce_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        //final result
        Utility.SetBuffer(Reduce_kernel, "reduceInput", hapticOutputBuffer);
        Utility.Dispatch(Reduce_kernel, Mathf.CeilToInt(population / (256f * 256f)), 1, 1);

        hapticOutputBuffer.GetData(result, 0, 0, 3);
        hapticInteractionPointDensity = result[0];
        //Debug.Log(hapticInteractionPointDensity);

        Simulation.SetFloat("hapticInteractionPointDensity", hapticInteractionPointDensity);
        Simulation.SetBuffer(Calculate_f_pressure_kernel, "particles", unsortedParticlesBuffer);
        Simulation.SetBuffer(Calculate_f_pressure_kernel, "hapticOutputForce", hapticOutputBuffer);
        Simulation.Dispatch(Calculate_f_pressure_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        //Simulation.SetVector("hapticInteractionPointDensity", hapticInVelocity);
        //Simulation.SetBuffer(Calculate_f_viscosity_kernel, "particles", unsortedParticlesBuffer);
        //Simulation.SetBuffer(Calculate_f_viscosity_kernel, "hapticOutputForce", hapticOutputBuffer);
        //Simulation.Dispatch(Calculate_f_viscosity_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        //per block result
        Utility.SetBuffer(Reduce_kernel, "reduceInput", hapticOutputBuffer);
        Utility.Dispatch(Reduce_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        //final result
        Utility.SetBuffer(Reduce_kernel, "reduceInput", hapticOutputBuffer);
        Utility.Dispatch(Reduce_kernel, Mathf.CeilToInt(population / (256f * 256f)), 1, 1);

        hapticOutputBuffer.GetData(result, 0, 0, 3);
        hapticPressureForce = new Vector3(result[0], result[1], result[2]);
        //Debug.Log(hapticPressureForce);


        Simulation.SetBuffer(Apply_f_pressure_kernel, "particles", unsortedParticlesBuffer);
        //reset the histogram buffer value to 0
        Simulation.SetBuffer(Apply_f_pressure_kernel, "globalHistogram", globalHistogramBuffer);
        Simulation.Dispatch(Apply_f_pressure_kernel, Mathf.CeilToInt(population / 256f), 1, 1);

        
        //indirect call to render
        instanceMaterial.SetBuffer("_Properties", meshPropertiesBuffer);
        Graphics.DrawMeshInstancedIndirect(sphereMesh, 0, instanceMaterial, bounds, argsBuffer);
    }
}
