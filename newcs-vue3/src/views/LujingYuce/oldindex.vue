<template>
    <div>
      <div class="header">合成路线智能预测</div>
  
      <div class="tab-container">
        <div class="tab">创建任务</div>
        <div class="tab">最近结果</div>
      </div>
  
      <div class="container">
        <div class="column">
          <div class="section upload-section">
            <h3>上传信息</h3>
            <input type="text" v-model="smiles" placeholder="输入SMILES字符串">
            <button class="submit-button" @click="submitTask">Start Explore</button>
          </div>
  
          <div class="section task-name-section">
            <h3>任务名称</h3>
            <input type="text" style="width: 100%; height: 30px;" v-model="taskName">
            <button class="submit-button" @click="submitTask">提交信息</button>
          </div>
        </div>
  
        <div class="column" v-if="predictionResults">
          <div class="section prediction-results">
            <h3>预测结果</h3>
            <p>目标分子: {{ predictionResults.target }}</p>
            <div>
              <h4>分子图像</h4>
              <!-- 图片路径需要带上/static -->
              <img :src="'http://localhost:5003/' + predictionResults.image_paths1[0]" alt="Target Image" />
            </div>
            <div v-for="(score, index) in predictionResults.scores_list" :key="index">
              <h4>可能的反应底物 {{ index + 1 }} - 概率: {{ score }}</h4>
              <div v-for="(path, imgIndex) in getImagePaths(index)" :key="imgIndex">
                <!-- 图片路径需要带上/static -->
                <img :src="'http://localhost:5003/' + path" :alt="'Result ' + (index + 1)" />
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  </template>
  
  <script>
  import axios from "axios";
  
  export default {
    data() {
      return {
        smiles: "",
        taskName: "",
        predictionResults: null,
      };
    },
    methods: {
      async submitTask() {
        try {
          const response = await axios.post("http://localhost:5003/predict", {
            smiles: this.smiles,
          });
  
          if (response.data.success) {
            this.predictionResults = response.data;
          } else {
            console.error("Prediction failed:", response.data.error);
          }
        } catch (error) {
          console.error("Error submitting task:", error);
        }
      },
      getImagePaths(index) {
        if (index === 0) {
          return this.predictionResults.image_paths1;
        } else if (index === 1) {
          return this.predictionResults.image_paths2;
        } else if (index === 2) {
          return this.predictionResults.image_paths3;
        }
        return [];
      },
    },
  };
  </script>
  
  <style>
  .header {
    background-color: #7abaff;
    padding: 10px;
    text-align: center;
    font-size: 1.5em;
    font-weight: bold;
  }
  .tab-container {
    display: flex;
    background-color: #e0e0e0;
    padding: 10px;
  }
  .tab {
    margin-right: 20px;
    padding: 10px;
    cursor: pointer;
  }
  .tab:hover {
    background-color: #d0d0d0;
  }
  .container {
    display: flex;
    flex-wrap: wrap;
    gap: 10px;
    padding: 20px;
  }
  .column {
    flex: 1;
    display: flex;
    flex-direction: column;
    gap: 10px;
  }
  .section {
    border: 1px solid #cccccc;
    padding: 10px;
    flex: 1;
  }
  .submit-button {
    background-color: #7abaff;
    border: none;
    padding: 10px;
    cursor: pointer;
    margin-top: 10px;
  }
  .submit-button:hover {
    background-color: #6aa2e0;
  }
  </style>
  