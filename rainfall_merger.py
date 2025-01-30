import numpy as np
import xarray as xr
import geopandas as gpd
import logging
from pykrige.ok import OrdinaryKriging
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

# Configure logging for error tracking
logging.basicConfig(level=logging.INFO)

def merge_gauge_satellite(gauge_points, satellite_grid, dem_grid):
    """
    Merges gauge data with satellite/DEM using Regression Kriging.
    
    Parameters:
    gauge_points (geopandas.GeoDataFrame): Must contain columns 'precip', 'lon', 'lat'.
    satellite_grid (xarray.DataArray): Satellite precipitation data (same CRS/resolution as DEM).
    dem_grid (xarray.DataArray): Digital Elevation Model data.
    
    Returns:
    xarray.DataArray: Merged precipitation field (mm/day) with the same grid as input rasters.
    
    Raises:
    ValueError: If input validation fails (missing columns, CRS mismatch, empty data).
    """
    try:
        # ----------------------
        # 1. Input Validation
        # ----------------------
        required_columns = ['precip', 'lon', 'lat']
        if not all(col in gauge_points.columns for col in required_columns):
            raise ValueError(f"gauge_points must contain columns: {required_columns}")
        if satellite_grid.rio.crs != dem_grid.rio.crs:
            raise ValueError("Satellite and DEM grids must share the same CRS")
        if gauge_points.empty or satellite_grid.size == 0 or dem_grid.size == 0:
            raise ValueError("Input data cannot be empty")

        # ----------------------
        # 2. Data Preparation
        # ----------------------
        lons = gauge_points['lon'].values
        lats = gauge_points['lat'].values
        precip = gauge_points['precip'].values

        # Extract satellite and DEM values at gauge locations
        sat_at_gauges = satellite_grid.sel(lon=lons, lat=lats, method="nearest").values
        dem_at_gauges = dem_grid.sel(lon=lons, lat=lats, method="nearest").values

        # Check for NaNs at gauge locations
        if np.isnan(sat_at_gauges).any() or np.isnan(dem_at_gauges).any():
            raise ValueError("NaN values detected at gauge locations in satellite/DEM data")

        # ----------------------
        # 3. Regression Model
        # ----------------------
        # Prepare features (satellite + DEM) and target (gauge precip)
        X = np.column_stack([sat_at_gauges, dem_at_gauges])
        y = precip

        # Split data into training/testing sets
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
        
        # Train linear regression model
        model = LinearRegression()
        model.fit(X_train, y_train)
        
        # Validate regression performance
        y_pred = model.predict(X_test)
        rmse = mean_squared_error(y_test, y_pred, squared=False)
        logging.info(f"Regression RMSE: {rmse:.2f} mm/day")

        # ----------------------
        # 4. Predict Trend
        # ----------------------
        # Flatten satellite and DEM grids for prediction
        sat_flat = satellite_grid.values.ravel()
        dem_flat = dem_grid.values.ravel()
        X_full = np.column_stack([sat_flat, dem_flat])
        
        # Predict trend (baseline precipitation) across the entire grid
        trend = model.predict(X_full).reshape(satellite_grid.shape)

        # ----------------------
        # 5. Krige Residuals
        # ----------------------
        # Calculate residuals (gauge - predicted trend)
        residuals = y - model.predict(X)
        
        # Configure Kriging (auto-fit variogram)
        ok = OrdinaryKriging(
            lons, lats, residuals,
            variogram_model='auto',  # Automatically fit variogram
            verbose=False  # Disable pykrige's internal prints
        )
        
        # Define output grid (same as satellite/DEM)
        grid_lon = satellite_grid.lon.values
        grid_lat = satellite_grid.lat.values
        
        # Interpolate residuals
        z_residual, _ = ok.execute('grid', grid_lon, grid_lat)

        # ----------------------
        # 6. Merge and Return
        # ----------------------
        # Combine trend and residuals
        merged = trend + z_residual.T  # Transpose to match (lat, lon) order
        
        # Wrap as xarray DataArray
        merged_da = xr.DataArray(
            merged,
            coords={"lat": satellite_grid.lat, "lon": satellite_grid.lon},
            dims=["lat", "lon"],
            attrs={"units": "mm/day", "method": "Regression Kriging"}
        )
        
        return merged_da

    except KeyError as e:
        logging.error(f"Missing required column: {e}")
        return None
    except ValueError as e:
        logging.error(f"Data validation error: {e}")
        return None
    except Exception as e:
        logging.error(f"Unexpected error during merging: {e}", exc_info=True)
        return None

# Example Usage
if __name__ == "__main__":
    # Load gauge data (point locations)
    gauge_points = gpd.read_file('gauges.shp')  # Requires 'precip', 'lon', 'lat'
    
    # Load satellite and DEM data (ensure same CRS/resolution)
    satellite = xr.open_dataarray('sat.tif').rio.write_crs("EPSG:4326")
    dem = xr.open_dataarray('dem.tif').rio.write_crs("EPSG:4326")
    
    # Run merging
    merged_precip = merge_gauge_satellite(gauge_points, satellite, dem)
    
    # Save output
    if merged_precip is not None:
        merged_precip.to_netcdf("merged_precipitation.nc")
        logging.info("Merged precipitation saved to merged_precipitation.nc")