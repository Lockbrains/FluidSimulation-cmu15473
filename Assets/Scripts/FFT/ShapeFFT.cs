using System;
using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using UnityEngine;

/// <summary>
/// For generating the Wave Shape using FFT.
/// </summary>
public class ShapeFFT : MonoBehaviour
{
    [SerializeField] private int resolution = 128;
    [SerializeField] private float frequency = 5f;

    [Header("For Spectrum Visualization")]
    [SerializeField] private MeshRenderer spectrumPlaneMR;
     
    private float[,] spectrum;
    private int spectrumSize;
    private Texture2D spectrumTexture;
    
    // Start is called before the first frame update
    void Start()
    {
        spectrum = new float[resolution, resolution];
        spectrumTexture = new Texture2D(resolution, resolution);
        //InitializeSpectrum();
    }

    // Update is called once per frame
    void Update()
    {
        InitializeSpectrum();
        GenerateSpectrum();
        PerformFFT();
    }
    
    private void InitializeSpectrum()
    {
        for (int x = 0; x < resolution; x++)
        {
            for (int y = 0; y < resolution; y++)
            {
                spectrum[x, y] = spectrum[x, y] = Mathf.Sin(frequency * x + Time.time) * Mathf.Cos(frequency * y + Time.time); 
            }
        }

        GenerateTextureFromSpectrum();
    }

    private void GenerateTextureFromSpectrum()
    {
        int size = spectrum.GetLength(0);

        if (!spectrumTexture)
        {
            spectrumTexture = new Texture2D(size, size);
        }
        
        Color[] pixels = new Color[size * size];
        
        for (int x = 0; x < size; x++)
        {
            for (int y = 0; y < size; y++)
            {
                float value = spectrum[x, y];
                pixels[x + y * size] = new Color(value, value, value, 1.0f);
            }
        }
        
        spectrumTexture.SetPixels(pixels);
        spectrumTexture.Apply();
    }

    public void GenerateSpectrum()
    {
        spectrumPlaneMR.material.mainTexture = spectrumTexture;
    }

    private void PerformFFT()
    {
        Complex[,] fftOutput = FFT2D(spectrum);
        
        int size = fftOutput.GetLength(0);

        for (int x = 0; x < size; x++)
        {
            for (int y = 0; y < size; y++)
            {
                float magnitude = (float) fftOutput[x, y].Magnitude;
                magnitude = Mathf.Clamp01(magnitude / 10f);
                
                spectrumTexture.SetPixel(x, y, new Color(magnitude, magnitude, magnitude, 1.0f));
            }
        }
        
        spectrumTexture.Apply();
    }
    
    #region 1D FFT
    public static Complex[] FFT1D(Complex[] input)
    {
        int N = input.Length;
        if (N <= 1) return input;

        var even = new Complex[N / 2];
        var odd  = new Complex[N / 2];
        for (int i = 0; i < N / 2; i++)
        {
            even[i] = input[i * 2];
            odd [i] = input[i * 2 + 1];
        }

        var fftEven = FFT1D(even);
        var fftOdd = FFT1D(odd);

        var result = new Complex[N];
        for (int k = 0; k < N / 2; k++)
        {
            Complex t = Complex.FromPolarCoordinates(1.0, -2.0 * Math.PI * k / N) * fftOdd[k];
            result[k] = fftEven[k] + t;
            result[k + N / 2] = fftEven[k] - t;
        }

        return result;
    }

    #endregion
    
    #region 2D FFT

    private Complex[,] FFT2D(float[,] input)
    {
        int size = input.GetLength(0);
        Complex[,] complexData = new Complex[size, size];

        // initialize input data
        for (int x = 0; x < size; x++)
        {
            for (int y = 0; y < size; y++)
            {
                complexData[x, y] = new Complex(input[x, y], 0);
            }
        }
        
        // run FFT1D for every row
        for (int y = 0; y < size; y++)
        {
            Complex[] row = new Complex[size];
            for (int x = 0; x < size; x++)
            {
                row[x] = complexData[x, y];
            }
            
            Complex[] fftRow = FFT1D(row);
            for (int x = 0; x < size; x++)
            {
                complexData[x, y] = fftRow[x];
            }
        }
        
        // run FFT1D for every col
        for (int x = 0; x < size; x++)
        {
            Complex[] column = new Complex[size];
            for (int y = 0; y < size; y++)
            {
                column[y] = complexData[x, y];
            }
            
            Complex[] fftColumn = FFT1D(column);
            for (int y = 0; y < size; y++)
            {
                complexData[x, y] = fftColumn[y]; 
            }
        }
        
        return complexData;
    }
    
    #endregion
} 
