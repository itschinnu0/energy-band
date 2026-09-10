"""Streamlit application entry point."""
import streamlit as st
import numpy as np
import scipy
import matplotlib.pyplot as plt

def main():
    st.title("Energy Band Diagram Simulator")
    st.subheader("Phase 1: Environment Spike")
    
    st.write("If you can see this, Streamlit is working!")
    st.write(f"NumPy version: {np.__version__}")
    st.write(f"SciPy version: {scipy.__version__}")
    
    fig, ax = plt.subplots()
    ax.plot([1, 2, 3], [1, 4, 9])
    st.pyplot(fig)

if __name__ == '__main__':
    main()
