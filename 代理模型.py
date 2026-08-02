import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

class SharedEncoder(nn.Module):
    def __init__(self, input_dim=17):
        super(SharedEncoder, self).__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(256, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(256, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.Dropout(0.3)
        )

    def forward(self, x):
        return self.net(x)

class BranchNetwork(nn.Module):
    def __init__(self, output_dim):
        super(BranchNetwork, self).__init__()
        self.net = nn.Sequential(
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Linear(32, output_dim)
        )

    def forward(self, x):
        return self.net(x)

class AircraftRunwaySurrogateModel(nn.Module):
    def __init__(self, input_dim=17):
        super(AircraftRunwaySurrogateModel, self).__init__()
        self.encoder = SharedEncoder(input_dim)
        # 总能耗分支
        self.total_energy_branch = BranchNetwork(output_dim=1)
        # 飞机分项 (缓冲器, 轮胎)
        self.aircraft_branch = BranchNetwork(output_dim=2)
        # 跑道分项 (面层, 基层, 道基)
        self.runway_branch = BranchNetwork(output_dim=3)

    def forward(self, x):
        features = self.encoder(x)
        return (self.total_energy_branch(features), 
                self.aircraft_branch(features), 
                self.runway_branch(features))

class PhysicsInformedLoss(nn.Module):
    def __init__(self, alpha_weights, lambda1=1.0, lambda2=1.0): # 权重取值可调整
        super(PhysicsInformedLoss, self).__init__()
        #  注册不需要梯度的权重参数
        self.register_buffer('alpha', alpha_weights)
        self.lambda1 = lambda1
        self.lambda2 = lambda2
        self.mse = nn.MSELoss()

    def forward(self, p_tot, p_air, p_run, t_tot, t_air, t_run):

        p_comps = torch.cat([p_air, p_run], dim=1)
        t_comps = torch.cat([t_air, t_run], dim=1)

        # 1. 分项能耗损失 
        huber_raw = F.huber_loss(p_comps, t_comps, reduction='none', delta=1.0)
        l_components = torch.mean(torch.sum(huber_raw * self.alpha, dim=1))

        # 2. 总能耗损失 (MSE) 
        l_total = self.mse(p_tot, t_tot)

        # 3. 物理守恒约束 
        sum_comps_pred = torch.sum(p_comps, dim=1, keepdim=True)
        l_physics = torch.mean(torch.abs(p_tot - sum_comps_pred))
       
        loss = l_components + self.lambda1 * l_total + self.lambda2 * l_physics
        return loss, l_components, l_total, l_physics


def load_and_process_data(csv_path, batch_size=64):

    df = pd.read_csv(csv_path)
    data = df.values.astype(np.float32) 

    X_raw = data[:, :17]   
    Y_total = data[:, 17:18]       
    Y_aircraft = data[:, 18:20]    
    Y_runway = data[:, 20:23]      
    
    Y_all = np.hstack([Y_total, Y_aircraft, Y_runway]) 

    # 划分数据集 (70% 训练, 15% 验证, 15% 测试) 
    X_train, X_temp, Y_train, Y_temp = train_test_split(
        X_raw, Y_all, test_size=0.30, random_state=42
    )

    X_val, X_test, Y_val, Y_test = train_test_split(
        X_temp, Y_temp, test_size=0.50, random_state=42
    )

    # 数据标准化
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    X_test = scaler.transform(X_test)

    # 计算动态权重
    train_components = Y_train[:, 1:] 
    means = np.mean(train_components, axis=0) + 1e-6
    alphas_numpy = 1.0 / means
    alphas_numpy = alphas_numpy / alphas_numpy.mean()
    
    print("-" * 30)
    print("自动计算的 Loss 权重 (Alphas):")
    print(f"缓冲器/轮胎: {alphas_numpy[0]:.4f}, {alphas_numpy[1]:.4f}")
    print(f"面层/基层/道基: {alphas_numpy[2]:.4f}, {alphas_numpy[3]:.4f}, {alphas_numpy[4]:.4f}")
    print("-" * 30)

    def create_loader(x, y):        
        t_tot = torch.tensor(y[:, 0:1])
        t_air = torch.tensor(y[:, 1:3])
        t_run = torch.tensor(y[:, 3:6])
        dataset = TensorDataset(torch.tensor(x), t_tot, t_air, t_run)
        return DataLoader(dataset, batch_size=batch_size, shuffle=(x is X_train))

    train_loader = create_loader(X_train, Y_train)
    val_loader = create_loader(X_val, Y_val)
    test_loader = create_loader(X_test, Y_test)

    return train_loader, val_loader, test_loader, alphas_numpy

def main():

    CSV_PATH = 'data.csv'  
    DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    LR = 3e-4 
    EPOCHS = 400
    
    try:
        train_loader, val_loader, test_loader, alpha_vals = load_and_process_data(CSV_PATH)
    except FileNotFoundError:
        print("错误: 找不到 data.csv，请确保文件在当前目录下。")
        return

    model = AircraftRunwaySurrogateModel(input_dim=17).to(DEVICE)
    
    alpha_tensor = torch.tensor(alpha_vals, dtype=torch.float32).to(DEVICE)
    criterion = PhysicsInformedLoss(alpha_weights=alpha_tensor).to(DEVICE)
    
    optimizer = optim.Adam(model.parameters(), lr=LR, weight_decay=1e-5) # [cite: 60]
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=EPOCHS)
    
    print("开始训练...")
    for epoch in range(EPOCHS):
        model.train()
        train_loss = 0.0
        
        for batch_x, batch_tot, batch_air, batch_run in train_loader:
            batch_x = batch_x.to(DEVICE)
            batch_tot = batch_tot.to(DEVICE)
            batch_air = batch_air.to(DEVICE)
            batch_run = batch_run.to(DEVICE)
            
            optimizer.zero_grad()
            
            p_tot, p_air, p_run = model(batch_x)
            
            loss, _, _, _ = criterion(p_tot, p_air, p_run, batch_tot, batch_air, batch_run)
            
            loss.backward()
            optimizer.step()
            train_loss += loss.item()
            
        scheduler.step()
        
        if epoch % 50 == 0:
            model.eval()
            val_phy_violation = 0.0
            with torch.no_grad():
                for bx, btot, bair, brun in val_loader:
                    bx, btot, bair, brun = bx.to(DEVICE), btot.to(DEVICE), bair.to(DEVICE), brun.to(DEVICE)
                    pt, pa, pr = model(bx)
                    _, _, _, l_phy = criterion(pt, pa, pr, btot, bair, brun)
                    val_phy_violation += l_phy.item()
            
            avg_phy = val_phy_violation / len(val_loader)
            print(f"Epoch {epoch} | Train Loss: {train_loss/len(train_loader):.4f} | Val Phy Violation: {avg_phy:.4f}")

    print("训练完成。")

if __name__ == "__main__":
    main()