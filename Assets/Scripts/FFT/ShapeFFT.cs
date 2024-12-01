using System;
using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using Unity.VisualScripting;
using UnityEngine;
using Vector2 = UnityEngine.Vector2;

/// <summary>
/// For generating the Wave Shape using FFT.
/// Creator: Lingheng Tony Tao
/// </summary>
public class ShapeFFT : MonoBehaviour
{
    [SerializeField] private int N = 256; // N
    float log2N => Mathf.Log(N, 2);
    [SerializeField] private float L = 1000f; // Patch Size
    [SerializeField] private float A = 0.01f; // Adjustment
    [SerializeField] private float V = 40f; // Wind Speed
    [SerializeField] private Vector2 windDirection = new Vector2(1, 1).normalized; 
    [SerializeField] private float frequency = 5f;

    [Header("For Spectrum Visualization")]
    [SerializeField] private MeshRenderer h0kRenderer;
    [SerializeField] private MeshRenderer h0MinusKRenderer;
    [SerializeField] private MeshRenderer butterflyTextureRenderer;
    [SerializeField] private MeshRenderer hktDyRenderer, hktDxRenderer, hktDzRenderer;
    [SerializeField] private MeshRenderer pingpongRenderer0, pingpongRenderer1;
    [SerializeField] private MeshRenderer pingpongRendererX0, pingpongRendererX1;
    [SerializeField] private MeshRenderer pingpongRendererZ0, pingpongRendererZ1;
    [SerializeField] private MeshRenderer foamTexture;
    [SerializeField] private MeshRenderer FFTResultRenderer;
    
    [Header("For Water Rendering")]
    [SerializeField] private MeshRenderer watermeshRenderer;
    [SerializeField] private float foamThreshold = 0.5f;
    [SerializeField] private float foamDecayRate = 0.1f;
    

    [Header("Compute Shaders")]
    [Tooltip("Don't use CPU. Will Crash.")]
    [SerializeField] private bool useGPU = true;
    [SerializeField] private ComputeShader fftInitShader;
    [SerializeField] private ComputeShader fourierComponentShader;
    [SerializeField] private ComputeShader butterflyTextureShader;
    [SerializeField] private ComputeShader butterflyComputeShader;
    [SerializeField] private ComputeShader inversionComputeShader;
    [SerializeField] private ComputeShader foamComputeShader;
    // Private Variables 
    private readonly float gravity = 9.81f;
    
    // For CPU
    private Complex[,] h0k;
    private Complex[,] h0MinusK;

    private Texture2D h0kTexture;
    private Texture2D h0MinusKTexture;
    private Texture2D butterflyTexture;
    
    private Texture2D pingpong0, pingpong1;
    
    // For GPU
    private int[] bitReversedIndices;
    
    private RenderTexture h0kRenderTexture;
    private RenderTexture h0MinusKRenderTexture;
    private RenderTexture butterflyRenderTexture;
    private RenderTexture pingpong0RT_Dx, pingpong1RT_Dx;
    private RenderTexture pingpong0RT_Dy, pingpong1RT_Dy;
    private RenderTexture pingpong0RT_Dz, pingpong1RT_Dz;
    private RenderTexture displacementRT;

    [SerializeField] private Texture2D perlinNoise;
    
    private RenderTexture hktDy, hktDx, hktDz;

    private RenderTexture foamRT;
    
    private int hktHandle;
    private int butterflyTextureHandle;
    private int butterflyComputeHandle;
    private int inversionHandle;
    private int foamHandle;
    
    private float time = 0f;

    void Start()
    {
        if (useGPU)
        {
            InitializeRenderTextures();
            InitializeSpectrumOnGPU();
            InitializeOtherComputeShaders();
        }
        else
        {
            InitializeSpectrum();

            int[] bitReversedIndices = GenerateBitReversedIndices(N);
            InitializeButterflyTexture(N, bitReversedIndices);

            InitializePingPongTextures(N);
        
            butterflyTextureRenderer.material.SetTexture("_MainTex", butterflyTexture);
            butterflyTexture.filterMode = FilterMode.Point;
        
            h0kTexture = GenerateTexture(h0k);
            h0kTexture.filterMode = FilterMode.Point;
        
            h0MinusKTexture = GenerateTexture(h0MinusK);
            h0MinusKTexture.filterMode = FilterMode.Point; 
        
            h0kRenderer.material.SetTexture("_MainTex", h0kTexture);
            h0MinusKRenderer.material.SetTexture("_MainTex", h0MinusKTexture);
        }
    }
    // Update is called once per frame
    void Update()
    {
        int currentPingpong = 0;
        time += Time.deltaTime;

        if (useGPU)
        {
            // Update Fourier Components
            fourierComponentShader.SetFloat("Time", time);

            int threadGroups = Mathf.CeilToInt(N / 16.0f);
            fourierComponentShader.Dispatch(hktHandle, threadGroups, threadGroups, 1);
            
            hktDyRenderer.material.SetTexture("_MainTex", hktDy);
            hktDzRenderer.material.SetTexture("_MainTex", hktDz);
            hktDxRenderer.material.SetTexture("_MainTex", hktDx);
            
            Graphics.Blit(hktDy, pingpong0RT_Dy);
            Graphics.Blit(hktDx, pingpong0RT_Dx);
            Graphics.Blit(hktDz, pingpong0RT_Dz);
            
            // horizontal fft
            for (int stage = 0; stage < log2N; stage++)
            {
                butterflyComputeShader.SetInt("stage", stage);
                butterflyComputeShader.SetInt("pingpong", currentPingpong);
                butterflyComputeShader.SetInt("direction", 0);
                
                DispatchButterflyComputeShader();
                
                currentPingpong = 1 - currentPingpong;
            }
            
            // vertical fft
            for (int stage = 0; stage < log2N; stage++)
            {
                butterflyComputeShader.SetInt("stage", stage);
                butterflyComputeShader.SetInt("pingpong", currentPingpong);
                butterflyComputeShader.SetInt("direction", 1);
                
                DispatchButterflyComputeShader();
                
                currentPingpong = 1 - currentPingpong;
            }
            
            pingpongRenderer0.material.SetTexture("_MainTex", pingpong0RT_Dy);
            pingpongRenderer1.material.SetTexture("_MainTex", pingpong1RT_Dy);
            pingpongRendererX0.material.SetTexture("_MainTex", pingpong0RT_Dx);
            pingpongRendererX1.material.SetTexture("_MainTex", pingpong1RT_Dx);
            pingpongRendererZ0.material.SetTexture("_MainTex", pingpong0RT_Dz);
            pingpongRendererZ1.material.SetTexture("_MainTex", pingpong1RT_Dz);
            
            // Inversion
            inversionComputeShader.Dispatch(inversionHandle, threadGroups, threadGroups, 1);
            FFTResultRenderer.material.SetTexture("_MainTex", displacementRT);
            
            watermeshRenderer.material.SetTexture("_DisplacementX", hktDx);
            watermeshRenderer.material.SetTexture("_DisplacementZ", hktDz);
            watermeshRenderer.material.SetTexture("_HeightMap", displacementRT);
            
            // Foam
            foamComputeShader.Dispatch(foamHandle, threadGroups, threadGroups, 1);
            watermeshRenderer.material.SetTexture("_FoamMap", foamRT);
        }
        else
        {
            // Update Fourier Components
            Texture2D hktTexture = UpdateFourierComponents(h0k, h0MinusK, time);

            CopyTextureData(hktTexture, pingpong1);

            hktDyRenderer.material.SetTexture("_MainTex", hktTexture); 

            
            CopyTextureData(hktTexture, pingpong0);
            
            currentPingpong = PerformHorizontalFFT();
            
            currentPingpong = PerformVerticalFFT(currentPingpong);

            Complex[,] result = ExtractDataFromTexture(currentPingpong == 0 ? pingpong0 : pingpong1);
            Complex[,] invertedResult = InverseAndPermuteCPU(result, N);
            
            Texture2D displacementTexture = GenerateTexture(invertedResult);
            FFTResultRenderer.material.SetTexture("_MainTex", displacementTexture);
            
            // Visualize 1D FFT Result
            pingpongRenderer0.material.SetTexture("_MainTex", pingpong0);
            pingpongRenderer1.material.SetTexture("_MainTex", pingpong1);
        }
    }
    
    #region Phase 1 Initialization
    // Phase 1: Produces the spectrum textures containing \Tilde{h}_0(\bold{k}) and \Tilde{h}_0^*(-\bold{k})
    // Page 24, 4.2.3
    private (float, float) GenerateGaussianRandom()
    {
        float u1 = UnityEngine.Random.value;
        float u2 = UnityEngine.Random.value;

        float theta = Mathf.Sin(2.0f * Mathf.PI * u2); 
        float z0 = Mathf.Sqrt(-2.0f * Mathf.Log(u1)) * Mathf.Cos(2.0f * Mathf.PI * u2);
        float z1 = Mathf.Sqrt(-2.0f * Mathf.Log(u1)) * Mathf.Sin(2.0f * Mathf.PI * u2);

        return (z0, z1);
    }
    
    // Generating Random Noise Texture for GPU 
    private Texture2D GenerateNoiseTexture(int N)
    {
        Texture2D noiseTexture = new Texture2D(N, N, TextureFormat.ARGB32, false);
        noiseTexture.filterMode = FilterMode.Point;
        noiseTexture.wrapMode = TextureWrapMode.Repeat;

        Color[] noisePixels = new Color[N * N];
        for (int i = 0; i < noisePixels.Length; i++)
        {
            float u1 = UnityEngine.Random.value;
            noisePixels[i] = new Color(u1, u1, u1, 0);
        }
        noiseTexture.SetPixels(noisePixels);
        noiseTexture.Apply();

        return noiseTexture;
    }
    
    /// <summary>
    /// Initialize Spectrum (h0(k) and h0(-k) on CPU
    /// </summary>
    private void InitializeSpectrum()
    {
        float windLength = V * V / gravity; // L = V^2/g;
        h0k = new Complex[N, N];
        h0MinusK = new Complex[N,N];
        
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                // Wave Vector \bold{k}
                // Equation (4.2)
                float kx = (2 * Mathf.PI * i - Mathf.PI * N ) / L;
                float kz = (2 * Mathf.PI * j - Mathf.PI * N ) / L;
                Vector2 k = new Vector2(kx, kz);
                
                // Clamp the size 
                float kLength = k.magnitude;
                if (kLength < 1e-6f) kLength = 1e-6f;

                // Phillips Spectrum
                // Equation (3.5)
                float kDotW = Vector2.Dot(k.normalized, windDirection);
                float kDotW2 = kDotW * kDotW;
                float phillips = A / Mathf.Pow(kLength, 4) * Mathf.Exp(-1 / Mathf.Pow(kLength * windLength, 2));
                phillips *= kDotW2;

                // Suppressing very small waves
                if (kLength < 1e-3f)
                    phillips = 0;

                // Generate Gaussian Random Numbers
                var (gaussianReal, gaussianImag) = GenerateGaussianRandom();

                // Initialize h0(k) and h0(-k)
                h0k[i, j] = new Complex(gaussianReal * Mathf.Sqrt(phillips / 2), gaussianImag * Mathf.Sqrt(phillips / 2));
                h0MinusK[i, j] = new Complex(h0k[i, j].Real, -h0k[i, j].Imaginary);
            }
        }
    }

    /// <summary>
    /// Initialize Spectrum (h0(k) and h0(-k) on GPU
    /// </summary>
    private void InitializeSpectrumOnGPU()
    {
        // Initialize h0(k) render texture
        h0kRenderTexture = new RenderTexture(N, N, 0, RenderTextureFormat.ARGBFloat);
        h0kRenderTexture.filterMode = FilterMode.Point; // no interpolation
        h0kRenderTexture.enableRandomWrite = true; // writable
        h0kRenderTexture.wrapMode = TextureWrapMode.Repeat;
        h0kRenderTexture.Create();
        
        // Initialize h0(-k) render texture
        h0MinusKRenderTexture = new RenderTexture(N, N, 0, RenderTextureFormat.ARGBFloat);
        h0MinusKRenderTexture.filterMode = FilterMode.Point; // no interpolation
        h0MinusKRenderTexture.enableRandomWrite = true; // writable
        h0MinusKRenderTexture.wrapMode = TextureWrapMode.Repeat;
        h0MinusKRenderTexture.Create();
        
        // find compute shader
        int kernel = fftInitShader.FindKernel("FFTInit");
        
        // set compute shader params
        fftInitShader.SetTexture(kernel, "H0k", h0kRenderTexture);
        fftInitShader.SetTexture(kernel, "H0MinusK", h0MinusKRenderTexture);
        
        // pass in other params
        fftInitShader.SetInt("N", N);
        fftInitShader.SetFloat("L", L);
        fftInitShader.SetFloat("A", A);
        fftInitShader.SetFloat("V", V);
        fftInitShader.SetFloats("windDirection", new float[] { windDirection.x, windDirection.y });
        
        // Dispatch the shader
        int threadGroupSize = Mathf.CeilToInt(N / 16.0f);
        fftInitShader.Dispatch(kernel, threadGroupSize, threadGroupSize, 1);
        
        // visualization
        h0kRenderer.material.SetTexture("_MainTex", h0kRenderTexture);
        h0MinusKRenderer.material.SetTexture("_MainTex", h0MinusKRenderTexture);
    }

    private RenderTexture CreateRenderTexture(int width, int height, bool interpolate = false, bool trillinear = false) {
        RenderTexture rt = new RenderTexture(width, height, 0, RenderTextureFormat.ARGBFloat);
        if (!interpolate)
        {
            rt.filterMode = FilterMode.Point;
        }
        else if (!trillinear)
        {
            rt.filterMode = FilterMode.Bilinear;
        }
        else
        {
            rt.filterMode = FilterMode.Trilinear;
        }
        rt.enableRandomWrite = true;
        rt.Create();
        return rt;
    }
    private void InitializeRenderTextures()
    {
        hktDy = CreateRenderTexture(N, N);
        hktDx = CreateRenderTexture(N, N);
        hktDz = CreateRenderTexture(N, N);
        
        int stages = Mathf.RoundToInt(Mathf.Log(N, 2)); 
        butterflyRenderTexture = new RenderTexture(stages, N,0, RenderTextureFormat.ARGBFloat);
        butterflyRenderTexture.filterMode = FilterMode.Point;
        butterflyRenderTexture.enableRandomWrite = true;
        butterflyRenderTexture.Create();
        
        pingpong0RT_Dy = CreateRenderTexture(N, N);
        pingpong1RT_Dy = CreateRenderTexture(N, N);
        pingpong0RT_Dx = CreateRenderTexture(N, N);
        pingpong1RT_Dx = CreateRenderTexture(N, N);
        pingpong0RT_Dz = CreateRenderTexture(N, N);
        pingpong1RT_Dz = CreateRenderTexture(N, N);
        
        displacementRT = CreateRenderTexture(N, N, true);
        foamRT = CreateRenderTexture(N, N, true, true);
        foamTexture.material.SetTexture("_MainTex", foamRT);
    }

    private void InitializeOtherComputeShaders()
    {
        hktHandle = fourierComponentShader.FindKernel("FourierComponent");
        
        // hkt 
        fourierComponentShader.SetTexture(hktHandle, "TildeHktDy", hktDy);
        fourierComponentShader.SetTexture(hktHandle, "TildeHktDx", hktDx);
        fourierComponentShader.SetTexture(hktHandle, "TildeHktDz", hktDz);
        
        fourierComponentShader.SetTexture(hktHandle, "Tilde_h0k", h0kRenderTexture);
        fourierComponentShader.SetTexture(hktHandle, "Tilde_h0MinusK", h0MinusKRenderTexture);
        
        fourierComponentShader.SetInt("N", N);
        fourierComponentShader.SetFloat("L", L);
        
        fourierComponentShader.SetFloat("Time", time);

        int threadGroups = Mathf.CeilToInt(N / 16.0f);
        fourierComponentShader.Dispatch(hktHandle, threadGroups, threadGroups, 1);
            
        hktDyRenderer.material.SetTexture("_MainTex", hktDy);
        hktDzRenderer.material.SetTexture("_MainTex", hktDz);
        hktDxRenderer.material.SetTexture("_MainTex", hktDx);
        
        Graphics.Blit(hktDy, pingpong0RT_Dy);
        Graphics.Blit(hktDx, pingpong0RT_Dx);
        Graphics.Blit(hktDz, pingpong0RT_Dz);
        
        // Butterfly Texture
        bitReversedIndices = GenerateBitReversedIndices(N);
        
        butterflyTextureHandle = butterflyTextureShader.FindKernel("ButterflyTexture");
        butterflyTextureShader.SetTexture(butterflyTextureHandle, "butterflyTexture", butterflyRenderTexture);
        butterflyTextureShader.SetInt("N", N);
        
        ComputeBuffer bitReversedBuffer = new ComputeBuffer(bitReversedIndices.Length, sizeof(int));
        bitReversedBuffer.SetData(bitReversedIndices);
        butterflyTextureShader.SetBuffer(butterflyTextureHandle, "bit_reversed", bitReversedBuffer);
        
        int threadGroupY = Mathf.CeilToInt(N / 16.0f);
        int threadGroupX = Mathf.CeilToInt(Mathf.Log(N, 2) / 16.0f);
        butterflyTextureShader.Dispatch(butterflyTextureHandle, threadGroupX, threadGroupY, 1);
        
        bitReversedBuffer.Release();
        
        butterflyTextureRenderer.material.SetTexture("_MainTex", butterflyRenderTexture);
        
        // Butterfly Compute 
        butterflyComputeHandle = butterflyComputeShader.FindKernel("Butterfly");
        butterflyComputeShader.SetTexture(butterflyComputeHandle, "butterflyTexture", butterflyRenderTexture);
        butterflyComputeShader.SetTexture(butterflyComputeHandle, "pingpong0_Dy", pingpong0RT_Dy);  
        butterflyComputeShader.SetTexture(butterflyComputeHandle, "pingpong1_Dy", pingpong1RT_Dy);  
        butterflyComputeShader.SetTexture(butterflyComputeHandle, "pingpong0_Dx", pingpong0RT_Dx);  
        butterflyComputeShader.SetTexture(butterflyComputeHandle, "pingpong1_Dx", pingpong1RT_Dx);
        butterflyComputeShader.SetTexture(butterflyComputeHandle, "pingpong0_Dz", pingpong0RT_Dz);  
        butterflyComputeShader.SetTexture(butterflyComputeHandle, "pingpong1_Dz", pingpong1RT_Dz);
        
        // Inversion
        inversionHandle = inversionComputeShader.FindKernel("IP");
        inversionComputeShader.SetInt("N", N);
        inversionComputeShader.SetTexture(inversionHandle, "pingpong0_Dy", pingpong0RT_Dy);
        inversionComputeShader.SetTexture(inversionHandle, "pingpong1_Dy", pingpong1RT_Dy);
        inversionComputeShader.SetTexture(inversionHandle, "pingpong0_Dx", pingpong0RT_Dx);
        inversionComputeShader.SetTexture(inversionHandle, "pingpong1_Dx", pingpong1RT_Dx);
        inversionComputeShader.SetTexture(inversionHandle, "pingpong0_Dz", pingpong0RT_Dz);
        inversionComputeShader.SetTexture(inversionHandle, "pingpong1_Dz", pingpong1RT_Dz);
        inversionComputeShader.SetTexture(inversionHandle, "displacement", displacementRT);
        
        // Foam
        foamHandle = foamComputeShader.FindKernel("GenerateFoam");
        foamComputeShader.SetTexture(foamHandle, "displacement", displacementRT);
        foamComputeShader.SetTexture(foamHandle, "foamTexture", foamRT);
        foamComputeShader.SetTexture(foamHandle, "perlinNoise", perlinNoise);
        foamComputeShader.SetFloat("threshold", foamThreshold);
        foamComputeShader.SetFloat("decayRate", foamDecayRate);
        foamComputeShader.SetFloat("deltaTime", Time.deltaTime);
        foamComputeShader.SetInt("N", N);
    }
    private int[] GenerateBitReversedIndices(int N)
    {
        int[] indices = new int[N];
        int log2N = Mathf.RoundToInt(Mathf.Log(N, 2));

        for (int i = 0; i < N; i++)
        {
            int reversed = 0;
            int index = i;

            for (int j = 0; j < log2N; j++)
            {
                reversed = (reversed << 1) | (index & 1);
                index >>= 1;
            }

            indices[i] = reversed;
        }
        return indices;
    }
    
    private void InitializeButterflyTexture(int N, int[] bitReversedIndices)
    {
        int stages = Mathf.RoundToInt(Mathf.Log(N, 2)); // Total Stage (log_2(N))
        butterflyTexture = new Texture2D(stages, N, TextureFormat.RGBAFloat, false);
        
        for (int y = 0; y < butterflyTexture.height; y++)
        {
            for (int x = 0; x < butterflyTexture.width; x++)
            {
                butterflyTexture.SetPixel(x, y, new Color(0, 0, 0, 0)); 
            }
        }
        butterflyTexture.Apply();

        for (int x = 0; x < stages; x++)
        {
            int butterflySpan = 1 << x; // current span

            for (int y = 0; y < N; y++)
            {
                // twiddle factor
                // Equation (4.6)
                float k = (y * (float)N / (1 << (x + 1))) % N;
                float twiddleReal = Mathf.Cos(2.0f * Mathf.PI * k / N);
                float twiddleImag = Mathf.Sin(2.0f * Mathf.PI * k / N);

                bool butterflyWing = ((y % (1 << (x + 1))) < butterflySpan); // wing mark

                // int flippedStage = stages - x - 1;
                // first stage, bit reversed indices
                if (x == 0)
                {
                    if (butterflyWing)
                    {
                        butterflyTexture.SetPixel(x, y, new Color(
                            twiddleReal,
                            twiddleImag,
                            bitReversedIndices[y],
                            bitReversedIndices[y + 1]
                        ));
                    }
                    else
                    {
                        butterflyTexture.SetPixel(x, y, new Color(
                            twiddleReal,
                            twiddleImag,
                            bitReversedIndices[y - 1],
                            bitReversedIndices[y]
                        ));
                    }
                }
                else
                {
                    if (butterflyWing)
                    {
                        butterflyTexture.SetPixel(x, y, new Color(
                            twiddleReal,
                            twiddleImag,
                            y,
                            y + butterflySpan
                        ));
                    }
                    else
                    {
                        butterflyTexture.SetPixel(x, y, new Color(
                            twiddleReal,
                            twiddleImag,
                            y - butterflySpan,
                            y
                        ));
                    }
                }
            }
        }
        butterflyTexture.Apply();
    }
    
    private Texture2D GenerateTexture(Complex[,] data)
    {
        Texture2D texture = new Texture2D(N, N, TextureFormat.RGFloat, false);

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                float real = (float)data[i, j].Real;
                float imag = (float)data[i, j].Imaginary;
                // The real part of h0(k) is stored in the red channel
                // while the imaginary part is stored in the green channel
                texture.SetPixel(i, j, new Color(real, imag, 0));
            }
        }

        texture.Apply();
        return texture;
    }
    
    private void InitializePingPongTextures(int N)
    {
        pingpong0 = new Texture2D(N, N, TextureFormat.RGBAFloat, false);
        pingpong0.filterMode = FilterMode.Point; 
        pingpong0.wrapMode = TextureWrapMode.Clamp;
        
        pingpong1 = new Texture2D(N, N, TextureFormat.RGBAFloat, false);
        pingpong1.filterMode = FilterMode.Point; 
        pingpong1.wrapMode = TextureWrapMode.Clamp; 
        
        Color[] blackPixels = new Color[N * N];
        for (int i = 0; i < blackPixels.Length; i++)
        {
            blackPixels[i] = new Color(0, 0, 0, 0); 
        }
        pingpong0.SetPixels(blackPixels);
        pingpong0.Apply();
        pingpong1.SetPixels(blackPixels);
        pingpong1.Apply();
    }
    
    #endregion
    
    #region Phase 2 Height Amplitude Fourier Components
    Texture2D UpdateFourierComponents(Complex[,] h0k, Complex[,] h0MinusK, float time)
    {
        Texture2D hktTexture = new Texture2D(N, N, TextureFormat.RGFloat, false);

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                // Wave Vector \bold{k}
                // Equation (4.2)
                float kx = (2 * Mathf.PI * i - Mathf.PI * N ) / L;
                float kz = (2 * Mathf.PI * j - Mathf.PI * N ) / L;
                Vector2 k = new Vector2(kx, kz);

                float kLength = k.magnitude;
                if (kLength < 1e-6f) kLength = 1e-6f;

                // frequency omega(k)
                // Equation (3.7)
                float omega = Mathf.Sqrt(9.8f * kLength);

                // Fourier Component
                // Equation (3.8)
                Complex expPositive = h0k[i, j] * Complex.Exp(new Complex(0, omega * time));
                Complex expNegative = h0MinusK[i, j] * Complex.Exp(new Complex(0, -omega * time));
                Complex hkt = expPositive + expNegative;

                // Write into Texture
                hktTexture.SetPixel(i, j, new Color((float)hkt.Real, (float)hkt.Imaginary, 0));
            }
        }

        hktTexture.Apply();
        return hktTexture;
    }
    
    #endregion
    
    #region Phase 3 Horitontal 1D FFT
    private void DispatchButterflyComputeShader() {
        int threadGroups = Mathf.CeilToInt(N / 16.0f);
        butterflyComputeShader.Dispatch(butterflyComputeHandle, threadGroups, threadGroups, 1);
    }
    
    private Complex[,] ExtractDataFromTexture(Texture2D texture)
    {
        int N = texture.width;
        Complex[,] data = new Complex[N, N];

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                Color color = texture.GetPixel(i, j);
                data[i, j] = new Complex(color.r, color.g); 
            }
        }

        return data;
    }

    Complex LoadComplexFromTexture(Texture2D texture, int x, int y)
    {
        Color pixel = texture.GetPixel(x, y);
        return new Complex(pixel.r, pixel.g); 
    }
    
    void StoreComplexToTexture(Texture2D texture, int x, int y, Color c)
    {
        texture.SetPixel(x, y, c);
        texture.Apply(); 
    }
    
    void CopyTextureData(Texture2D source, Texture2D target)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                Color color = source.GetPixel(i, j);
                target.SetPixel(i, j, color);
            }
        }
        target.Apply();
    }
    
    private void HorizontalButterflies(int stage, int pingpong)
    {
        Color[] butterflyData = butterflyTexture.GetPixels(); 
        Complex[,] input = pingpong == 0 ? ExtractDataFromTexture(pingpong1) : ExtractDataFromTexture(pingpong0);
        Complex[,] output = new Complex[N, N];
        
        Color[] resultPixels = new Color[N * N];
        // Phase 3 Horitonzal 1D FFT
        for (int y = 0; y < N; y++)
        {
            for (int x = 0; x < N; x++)
            {
                Color data = butterflyData[stage * N + x];
                
                int px = Mathf.FloorToInt(data.b);
                int qx = Mathf.FloorToInt(data.a);
                
                Complex p, q;
                p = input[px, y];
                q = input[qx, y];
                Complex w = new Complex(data.r, data.g);

                Complex H = p + w * q;
                output[x, y] = H;
                resultPixels[y * N + x] = new Color((float)H.Real, (float)H.Imaginary, 0, 1);
            }
        }
        
        if (pingpong == 0)
        {
            pingpong1.SetPixels(resultPixels);
            pingpong1.Apply(); 
        }
        else
        {
            pingpong0.SetPixels(resultPixels);
            pingpong0.Apply(); 
        }
    }
    
    int PerformHorizontalFFT()
    {
        int stages = Mathf.RoundToInt(Mathf.Log(N, 2));
        int pingpong = 0;

        for (int stage = 0; stage < stages; stage++)
        {
            HorizontalButterflies(stage, pingpong);
            pingpong = 1 - pingpong;
        }

        return pingpong;
    }
    
    
    
    #endregion
    
    #region Phase 4 Vertical 1D FFT
    int PerformVerticalFFT(int initPingpong)
    {
        int stages = Mathf.RoundToInt(Mathf.Log(N, 2));
        int pingpong = initPingpong;

        for (int stage = 0; stage < stages; stage++)
        {
            VerticalButterflies(stage, pingpong);
            pingpong = 1 - pingpong; 
        }

        return pingpong;
    }
    
    private void VerticalButterflies(int stage, int pingpong)
    {
        Color[] butterflyData = butterflyTexture.GetPixels();
        Complex[,] input = pingpong == 0 ? ExtractDataFromTexture(pingpong1) : ExtractDataFromTexture(pingpong0);
        Complex[,] output = new Complex[N, N];

        Color[] resultPixels = new Color[N * N]; 
        
        for (int x = 0; x < N; x++)
        {
            for (int y = 0; y < N; y++)
            {
                Color data = butterflyData[stage * N + y];
                int py = Mathf.FloorToInt(data.b);
                int qy = Mathf.FloorToInt(data.a);
                
                Complex p = input[x, py];
                Complex q = input[x, qy];
                Complex w = new Complex(data.r, data.g);

                Complex H = p + w * q;
                
                output[x, y] = H;
                resultPixels[y * N + x] = new Color((float)H.Real, (float)H.Imaginary, 0, 1);
            }
        }
        
        if (pingpong == 0)
        {
            pingpong1.SetPixels(resultPixels);
            pingpong1.Apply(); 
        }
        else
        {
            pingpong0.SetPixels(resultPixels);
            pingpong0.Apply(); 
        }
    }
    
    #endregion
    
    #region Phase 5 Inversion
    private Complex[,] InverseAndPermuteCPU(Complex[,] input, int N)
    {
        Complex[,] output = new Complex[N, N];

        // 1. Define perms array for alternation
        float[] perms = { 1.0f, -1.0f };

        // 2. Apply scaling factor and permutation
        for (int y = 0; y < N; y++)
        {
            for (int x = 0; x < N; x++)
            {
                // Determine perm based on (x + y) % 2
                int index = (x + y) % 2;
                float perm = perms[index];

                // Apply scaling factor (1 / N^2) and permutation
                Complex h = input[x, y] * (1.0f / (N * N));
                output[x, y] = new Complex(h.Real * perm, h.Imaginary * perm);
            }
        }

        return output;
    }
    #endregion
} 
